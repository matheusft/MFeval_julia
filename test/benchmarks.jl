# ==============================================================================
# benchmarks.jl
#
# Phase 5 — Performance
#
# Five sections covering all spec requirements:
#   1. Scalar path        — target < 1 µs per evaluation
#   2. Allocation audit   — zero heap allocs on scalar hot path
#   3. Type stability     — zero Any / Union{} on full hot path
#   4. Batch throughput   — ns/eval for N = 100 / 1_000 / 10_000
#   5. Thread scaling     — parallel speedup vs ideal serial baseline
#
# Run standalone:
#   julia --project=. test/benchmarks.jl
#   julia --project=. --threads=4 test/benchmarks.jl
#
# NOT included in runtests.jl — benchmarks are not a correctness gate.
# ==============================================================================

include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
using .MFeval
using BenchmarkTools
using Printf
using InteractiveUtils

# ── Fixtures ──────────────────────────────────────────────────────────────────
const FIX = joinpath(@__DIR__, "fixtures")

const _bp61  = read_tir(joinpath(FIX, "MagicFormula61_Parameters.tir"))
const _bp62  = read_tir(joinpath(FIX, "MagicFormula62_Parameters.tir"))
const _bp52  = read_tir(joinpath(FIX, "MagicFormula52_Parameters.tir"))
const _bM111 = MFModes(111)

const _bNOM61 = MFInputs(_bp61.fnomin, -0.05, 0.1, 0.02, 0.0, _bp61.longvl)
const _bNOM62 = MFInputs(_bp62.fnomin, -0.05, 0.1, 0.02, 0.0, _bp62.longvl)
const _bNOM52 = MFInputs(_bp52.fnomin, -0.05, 0.1, 0.02, 0.0, _bp52.longvl)

const TARGET_NS = 1_000.0   # 1 µs

# ── Helpers ───────────────────────────────────────────────────────────────────
pass(s)  = "\e[32m PASS \e[0m  $s"
fail(s)  = "\e[31m FAIL \e[0m  $s"
warn(s)  = "\e[33m WARN \e[0m  $s"
head(s)  = "\e[1m$s\e[0m"

function hline(); println("─" ^ 72); end
function section(s); hline(); println(head(s)); hline(); end

# ── 1. SCALAR BENCHMARKS ──────────────────────────────────────────────────────

section("1 · Scalar path  (target < 1 µs)")

for (label, tp, inp) in (
        ("MF6.1  combined slip + camber", _bp61, _bNOM61),
        ("MF6.2  (iterative Rl solver)",  _bp62, _bNOM62),
        ("MF5.2",                          _bp52, _bNOM52),
        ("MF6.1  turn-slip (M112)",        _bp61,
            MFInputs(_bp61.fnomin, 0.0, 0.1, 0.0, 0.02, _bp61.longvl)),
    )
    b   = @benchmark mfeval($tp, $inp, $_bM111) evals=500 samples=200 seconds=3
    med = median(b).time
    mn  = minimum(b).time
    mx  = maximum(b).time
    msg = @sprintf "%-42s  %7.0f ns  (min %5.0f, max %6.0f)" label med mn mx
    println(med < TARGET_NS ? pass(msg) : (med < 2*TARGET_NS ? warn(msg) : fail(msg)))
end

# ── 2. ALLOCATION AUDIT ───────────────────────────────────────────────────────

section("2 · Allocation audit  (scalar path must be zero-alloc)")

for (label, tp, inp) in (
        ("MF6.1 scalar", _bp61, _bNOM61),
        ("MF6.2 scalar", _bp62, _bNOM62),
        ("MF5.2 scalar", _bp52, _bNOM52),
    )
    b      = @benchmark mfeval($tp, $inp, $_bM111) evals=500 samples=50
    allocs = b.allocs
    msg    = @sprintf "%-20s  %d allocation(s)" label allocs
    println(allocs == 0 ? pass(msg) : fail(msg))
end

# ── 3. TYPE STABILITY ─────────────────────────────────────────────────────────

section("3 · Type stability  (@code_warntype audit)")

function check_warntype(label, f, args...)
    buf = IOBuffer()
    code_warntype(IOContext(buf, :color => false), f, map(typeof, args))
    output = String(take!(buf))
    # Skip the signature header lines — only scan the body
    body = join(split(output, "\n")[3:end], "\n")
    has_any   = occursin("::Any",    body)
    has_union = occursin("::Union{", body)
    if has_any || has_union
        println(fail(@sprintf "%-42s  instability detected" label))
        for line in split(body, "\n")
            (occursin("::Any", line) || occursin("::Union{", line)) &&
                println("      → ", strip(line))
        end
    else
        println(pass(@sprintf "%-42s" label))
    end
end

check_warntype("mfeval scalar  MF6.1", mfeval, _bp61, _bNOM61, _bM111)
check_warntype("mfeval scalar  MF6.2", mfeval, _bp62, _bNOM62, _bM111)
check_warntype("mfeval scalar  MF5.2", mfeval, _bp52, _bNOM52, _bM111)

# ── 4. BATCH THROUGHPUT ───────────────────────────────────────────────────────

section("4 · Batch throughput  (ns/eval, N = 100 / 1_000 / 10_000)")

let
    Fnom = _bp61.fnomin
    Vnom = _bp61.longvl
    thr  = Threads.nthreads()

    for N in (100, 1_000, 10_000)
        mat = Matrix{Float64}(undef, N, 6)
        for i in 1:N
            mat[i, :] = [Fnom, 0.0, (i/N - 0.5)*0.6, 0.0, 0.0, Vnom]
        end
        out = Matrix{Float64}(undef, N, 30)

        b      = @benchmark mfeval!($out, $_bp61, $mat, $_bM111) evals=20 samples=100
        med_ns = median(b).time
        ns_per = med_ns / N
        msg    = @sprintf "N=%6d  %2d thread(s)   %7.0f ns/eval   total %6.2f ms" N thr ns_per (med_ns/1e6)
        println(ns_per < TARGET_NS ? pass(msg) : warn(msg))
    end
end

# ── 5. THREAD SCALING ─────────────────────────────────────────────────────────

section("5 · Thread scaling  (parallel speedup vs ideal serial)")

let
    nthreads = Threads.nthreads()
    if nthreads < 2
        println("  Skipped — start Julia with --threads=N (N ≥ 2) to measure scaling.")
    else
        N    = 10_000
        Fnom = _bp61.fnomin
        Vnom = _bp61.longvl
        mat  = Matrix{Float64}(undef, N, 6)
        for i in 1:N
            mat[i, :] = [Fnom, 0.0, (i/N - 0.5)*0.6, 0.0, 0.0, Vnom]
        end
        out = Matrix{Float64}(undef, N, 30)

        b_batch  = @benchmark mfeval!($out, $_bp61, $mat, $_bM111) evals=20 samples=100
        b_scalar = @benchmark mfeval($_bp61, $_bNOM61, $_bM111) evals=500 samples=100

        t_batch  = median(b_batch).time
        t_scalar = median(b_scalar).time
        ideal    = t_scalar * N
        speedup  = ideal / t_batch
        eff      = speedup / nthreads

        println("  Threads:           $nthreads")
        println("  Scalar median:     $(round(t_scalar, digits=0)) ns")
        println("  Batch median:      $(@sprintf "%.2f" t_batch/1e6) ms  (N=$N)")
        println("  Ideal serial est:  $(@sprintf "%.2f" ideal/1e6) ms")
        println("  Parallel speedup:  $(@sprintf "%.2f" speedup)×")
        println("  Thread efficiency: $(@sprintf "%.0f" eff*100)%")
        msg = @sprintf "efficiency %.0f%% (≥ 50%% target)" (eff*100)
        println(eff >= 0.5 ? pass(msg) : warn(msg))
    end
end

hline()
println("Done.")