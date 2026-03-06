# crossval_julia.jl
# ==========================================================================
# Cross-validation & benchmark suite for MFeval.jl
#
# Completely independent of MATLAB at runtime. Reads the CSV files that
# crossval_matlab.m produced (in test/crossval_data/) and compares every
# output column against Julia's own mfeval evaluation of the same inputs.
#
# Each CSV has the layout:
#   input columns (6 or 7), then 30 output columns
# The script detects which format by reading the header.
#
# Usage:
#   julia --project=. --threads=auto test/crossval_julia.jl
#
# Output:
#   Console: pass/fail per scenario, column-level error stats
#   test/crossval_perf_julia.csv    — Julia benchmark results
#   test/crossval_comparison.csv    — side-by-side table (if MATLAB perf CSV exists)
#   test/crossval_error_report.csv  — per-test error summary
# ==========================================================================

include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
using .MFeval
using BenchmarkTools
using Printf
using Statistics

# ── Configuration ─────────────────────────────────────────────────────────────
const FIX_DIR  = joinpath(@__DIR__, "fixtures")
const DATA_DIR = joinpath(@__DIR__, "crossval_data")
const PERF_CSV = joinpath(@__DIR__, "crossval_perf_julia.csv")
const COMP_CSV = joinpath(@__DIR__, "crossval_comparison.csv")

const TOL_DEFAULT  = 1e-6
const TOL_INST_KYA = 1e-3   # col 29 — finite-difference quantity

const COL_NAMES = ["Fx","Fy","Fz_out","Mx","My","Mz",
                   "kappa_out","alpha_out","gamma_out","phit_out","Vx_out","pressure_out",
                   "Re","rho","two_a","t","mux","muy","omega","Rl","two_b","Mzr",
                   "Cx","Cy","Cz","Kya","sigmax","sigmay","inst_Kya","Kxk"]

# ── Load TIR files ────────────────────────────────────────────────────────────
const p61 = read_tir(joinpath(FIX_DIR, "MagicFormula61_Parameters.tir"))
const p62 = read_tir(joinpath(FIX_DIR, "MagicFormula62_Parameters.tir"))
const p52 = read_tir(joinpath(FIX_DIR, "MagicFormula52_Parameters.tir"))
const M111 = MFModes(111)

println("=============================================================")
println("  Julia MFeval.jl Cross-Validation & Benchmark Suite")
println("  Julia $(VERSION), $(Threads.nthreads()) thread(s)")
println("=============================================================\n")

# ======================================================================
# CSV helpers
# ======================================================================

function load_csv_matrix(path::String)
    lines = readlines(path)
    header = String.(split(lines[1], ','))
    N = length(lines) - 1
    ncols = length(header)
    mat = Matrix{Float64}(undef, N, ncols)
    for (i, line) in enumerate(lines[2:end])
        mat[i, :] = parse.(Float64, split(line, ','))
    end
    return mat, header
end

"""
    load_crossval_csv(name) → (inputs, matlab_ref, n_input_cols) or nothing

Load a combined CSV from crossval_data/. Detects 6 or 7 input columns
by checking the header for a 'pressure' column. Returns the input matrix,
the 30-column MATLAB reference output matrix, and the number of input cols.
"""
function load_crossval_csv(name::String)
    path = joinpath(DATA_DIR, "$(name).csv")
    isfile(path) || return nothing
    mat, header = load_csv_matrix(path)
    # Detect input column count from header
    n_in = findfirst(==("Fx"), header) - 1  # "Fx" is the first output column
    if n_in === nothing || n_in < 6
        # Fallback: try "Fz_out" as the first output marker
        n_in_alt = findfirst(==("Fz_out"), header)
        if n_in_alt !== nothing
            # The output columns start at "Fx" which is 2 before "Fz_out"
            n_in = n_in_alt - 3  # Fx, Fy, then Fz_out
        end
        # Ultimate fallback: check if total cols = 36 (6+30) or 37 (7+30)
        if n_in === nothing || n_in < 6
            total_cols = size(mat, 2)
            n_in = total_cols - 30
        end
    end
    inputs     = mat[:, 1:n_in]
    matlab_ref = mat[:, n_in+1:n_in+30]
    return (inputs=inputs, ref=matlab_ref, n_in=n_in)
end

# ======================================================================
# Comparison engine
# ======================================================================

function compare_outputs(julia_out::Matrix{Float64}, matlab_ref::Matrix{Float64},
                         name::String; verbose::Bool=true)
    N = size(matlab_ref, 1)
    @assert size(julia_out) == size(matlab_ref) "Size mismatch: Julia=$(size(julia_out)), MATLAB=$(size(matlab_ref))"

    col_max_err  = zeros(30)
    col_mean_err = zeros(30)
    col_pass     = trues(30)
    n_fail       = 0

    for c in 1:30
        tol = (c == 29) ? TOL_INST_KYA : TOL_DEFAULT

        diffs = abs.(julia_out[:, c] .- matlab_ref[:, c])
        scales = max.(abs.(matlab_ref[:, c]), 1.0)
        rel_diffs = diffs ./ scales

        # Filter out NaN/Inf for statistics (they are handled separately)
        valid_diffs = filter(!isnan, diffs)
        valid_diffs = filter(!isinf, valid_diffs)
        
        col_max_err[c]  = isempty(valid_diffs) ? NaN : maximum(valid_diffs)
        col_mean_err[c] = isempty(valid_diffs) ? NaN : mean(valid_diffs)

        for i in 1:N
            jval = julia_out[i, c]
            mval = matlab_ref[i, c]
            
            # Handle NaN/Inf matching: both must be same special value
            if isnan(jval) && isnan(mval)
                continue  # NaN == NaN is acceptable
            elseif isinf(jval) && isinf(mval) && sign(jval) == sign(mval)
                continue  # +Inf == +Inf or -Inf == -Inf is acceptable
            elseif isnan(jval) || isnan(mval) || isinf(jval) || isinf(mval)
                col_pass[c] = false  # One is NaN/Inf, other isn't
                n_fail += 1
                continue
            end
            
            # Normal numeric comparison
            ref_val = abs(mval)
            err = (ref_val < 1e-10) ? diffs[i] : rel_diffs[i]
            if err > tol
                col_pass[c] = false
                n_fail += 1
            end
        end
    end

    all_pass = all(col_pass)

    if verbose
        status = all_pass ? "PASS ✓" : "FAIL ✗"
        @printf("  %-30s  %s  (%d rows)\n", name, status, N)
        if !all_pass
            for c in 1:30
                if !col_pass[c]
                    @printf("    col %2d %-14s  max_err=%.3e  mean_err=%.3e\n",
                            c, COL_NAMES[c], col_max_err[c], col_mean_err[c])
                end
            end
        end
    end

    return (pass=all_pass, N=N, col_max_err=col_max_err,
            col_mean_err=col_mean_err, col_pass=col_pass, n_fail=n_fail)
end

# ======================================================================
# Param lookup
# ======================================================================

function get_params(name::String)
    startswith(name, "mf61") && return p61
    startswith(name, "mf62") && return p62
    startswith(name, "mf52") && return p52
    error("Cannot determine MF version from name: $name")
end

# ======================================================================
# PART A: Cross-Validation
# ======================================================================
println("============ PART A: Cross-Validation ============\n")

# Check if crossval_data directory exists
if !isdir(DATA_DIR)
    println("  ⚠  Directory not found: $DATA_DIR")
    println("     Run crossval_matlab.m in MATLAB first, then copy")
    println("     the test/crossval_data/ folder into this project.\n")
end

# Standard test cases (all useMode = 111)
const STANDARD_TESTS = [
    "mf61_pure_lat", "mf62_pure_lat", "mf52_pure_lat",
    "mf61_pure_lon", "mf62_pure_lon", "mf52_pure_lon",
    "mf61_combined", "mf62_combined", "mf52_combined",
    "mf61_camber",   "mf62_camber",   "mf52_camber",
    "mf61_load_sweep",
    "mf61_pressure",
    "mf61_carpet_lat", "mf61_carpet_lon", "mf61_carpet_combined",
    "mf61_edge_cases",
    "mf61_velocity",
    "mf61_canonical", "mf62_canonical", "mf52_canonical",
]

# Wrap Part A in a function to avoid soft-scope issues
function run_validation()
    total_pass  = 0
    total_fail  = 0
    total_tests = 0
    total_rows  = 0
    results_summary = NamedTuple{(:name,:pass,:N,:max_err), Tuple{String,Bool,Int,Float64}}[]

    println("--- Standard tests (useMode = 111) ---\n")

    for name in STANDARD_TESTS
        data = load_crossval_csv(name)
        if data === nothing
            @printf("  %-30s  SKIP (no reference CSV)\n", name)
            continue
        end

        total_tests += 1
        p = get_params(name)
        julia_out = mfeval(p, data.inputs, M111)
        r = compare_outputs(julia_out, data.ref, name)
        total_rows += r.N

        if r.pass
            total_pass += 1
        else
            total_fail += 1
        end
        push!(results_summary, (name=name, pass=r.pass, N=r.N,
                                max_err=maximum(r.col_max_err)))
    end

    # ── useMode variation test ────────────────────────────────────────────
    println("\n--- useMode variation tests ---\n")

    mode_path = joinpath(DATA_DIR, "mf61_usemodes.csv")
    if isfile(mode_path)
        mat, header = load_csv_matrix(mode_path)
        # Layout: Fz,kappa,alpha,gamma,phit,Vx,useMode, then 30 output cols
        useModes    = Int.(mat[:, 7])
        matlab_ref  = mat[:, 8:37]
        canonical   = [4000.0, -0.05, 0.1, 0.02, 0.0, 16.7]
        nModes      = length(useModes)
        julia_out   = Matrix{Float64}(undef, nModes, 30)

        for (i, um) in enumerate(useModes)
            modes_i = MFModes(um)
            inp_i   = MFInputs(canonical[1], canonical[2], canonical[3],
                                canonical[4], canonical[5], canonical[6])
            julia_out[i, :] = to_vector(mfeval(p61, inp_i, modes_i))
        end

        total_tests += 1
        r = compare_outputs(julia_out, matlab_ref, "mf61_usemodes")
        total_rows += r.N
        if r.pass
            total_pass += 1
        else
            total_fail += 1
        end
        push!(results_summary, (name="mf61_usemodes", pass=r.pass, N=r.N,
                                max_err=maximum(r.col_max_err)))
    else
        @printf("  %-30s  SKIP (no reference CSV)\n", "mf61_usemodes")
    end

    # ── Summary ───────────────────────────────────────────────────────────
    println("\n--- Validation Summary ---")
    @printf("  Tests run:    %d\n", total_tests)
    @printf("  Tests passed: %d\n", total_pass)
    @printf("  Tests failed: %d\n", total_fail)
    @printf("  Total rows:   %d\n", total_rows)

    if total_fail > 0
        println("\n  Failed tests:")
        for r in results_summary
            if !r.pass
                @printf("    %-30s  max_err=%.3e  (%d rows)\n", r.name, r.max_err, r.N)
            end
        end
    end

    # ── Write error report ────────────────────────────────────────────────
    err_path = joinpath(@__DIR__, "crossval_error_report.csv")
    open(err_path, "w") do io
        println(io, "test_name,N,pass,max_error")
        for r in results_summary
            @printf(io, "%s,%d,%s,%.6e\n", r.name, r.N, r.pass ? "PASS" : "FAIL", r.max_err)
        end
    end
    println("  Error report: $err_path")

    return (total_pass=total_pass, total_fail=total_fail, total_tests=total_tests)
end

val_result = run_validation()


# ======================================================================
# PART B: Performance Benchmarks
# ======================================================================
println("\n============ PART B: Performance Benchmarks ============\n")

function run_benchmarks()
    perf_results = Tuple{String,Int,Float64,Float64,Float64,Float64,Float64}[]

    # ── B1. Scalar — all versions ─────────────────────────────────────────
    println("--- B1. Scalar evaluation (single point) ---")
    for (label, tp, inp) in [
        ("MF 5.2", p52, MFInputs(p52.fnomin, -0.05, 0.1, 0.02, 0.0, p52.longvl)),
        ("MF 6.1", p61, MFInputs(p61.fnomin, -0.05, 0.1, 0.02, 0.0, p61.longvl)),
        ("MF 6.2", p62, MFInputs(p62.fnomin, -0.05, 0.1, 0.02, 0.0, p62.longvl)),
    ]
        b = @benchmark mfeval($tp, $inp, $M111) evals=500 samples=500 seconds=5
        med = median(b).time
        mn  = minimum(b).time
        mx  = maximum(b).time
        p05 = quantile(b.times, 0.05)
        p95 = quantile(b.times, 0.95)
        @printf("  %-8s  median: %7.1f ns   min: %7.1f ns   max: %7.1f ns   p5: %7.1f   p95: %7.1f\n",
                label, med, mn, mx, p05, p95)
        push!(perf_results, ("scalar_$label", 1, med, mn, mx, p05, p95))
    end
    println()

    # ── B2. Batch — MF 6.1 ───────────────────────────────────────────────
    println("--- B2. Batch evaluation (MF 6.1, $(Threads.nthreads()) threads) ---")
    for Nb in [10, 50, 100, 500, 1000, 5000, 10000]
        mat = Matrix{Float64}(undef, Nb, 6)
        for i in 1:Nb
            mat[i, :] = [4000.0, 0.0, -0.3 + 0.6*(i-1)/max(Nb-1,1), 0.0, 0.0, 16.7]
        end
        out = Matrix{Float64}(undef, Nb, 30)
        b = @benchmark mfeval!($out, $p61, $mat, $M111) evals=20 samples=200 seconds=5
        med = median(b).time
        ns_per = med / Nb
        @printf("  N=%6d   median: %8.3f ms   (%7.1f ns/eval)\n", Nb, med/1e6, ns_per)
        push!(perf_results, ("batch_N$Nb", Nb, ns_per,
              minimum(b).time/Nb, maximum(b).time/Nb,
              quantile(b.times, 0.05)/Nb, quantile(b.times, 0.95)/Nb))
    end
    println()

    # ── B3. Pure lateral sweep N=200 ──────────────────────────────────────
    println("--- B3. Pure lateral sweep (N=200, Fz=3000 N) ---")
    nP = 200
    mat200 = Matrix{Float64}(undef, nP, 6)
    for i in 1:nP
        mat200[i, :] = [3000.0, 0.0, -0.3 + 0.6*(i-1)/(nP-1), 0.0, 0.0, 16.0]
    end
    out200 = Matrix{Float64}(undef, nP, 30)
    b200 = @benchmark mfeval!($out200, $p61, $mat200, $M111) evals=50 samples=1000 seconds=10
    med200 = median(b200).time
    @printf("  median: %.4f ms   min: %.4f ms   max: %.4f ms   (%.1f ns/eval)\n",
            med200/1e6, minimum(b200).time/1e6, maximum(b200).time/1e6, med200/nP)
    push!(perf_results, ("lat_sweep_200", 200, med200/nP,
          minimum(b200).time/nP, maximum(b200).time/nP,
          quantile(b200.times, 0.05)/nP, quantile(b200.times, 0.95)/nP))
    println()

    # ── B4. Multi-load carpet 4x200=800 ──────────────────────────────────
    println("--- B4. Multi-load lateral carpet (4 x 200 = 800) ---")
    Fz_vals = [2000.0, 4000.0, 6000.0, 8000.0]
    nAlpha = 200; nTotal = length(Fz_vals) * nAlpha
    mat_c = Matrix{Float64}(undef, nTotal, 6)
    row = 1
    for Fz in Fz_vals, i in 1:nAlpha
        mat_c[row, :] = [Fz, 0.0, -0.5 + 1.0*(i-1)/(nAlpha-1), 0.0, 0.0, 16.7]
        row += 1
    end
    out_c = Matrix{Float64}(undef, nTotal, 30)
    b_c = @benchmark mfeval!($out_c, $p61, $mat_c, $M111) evals=20 samples=200 seconds=5
    med_c = median(b_c).time
    @printf("  800 evals: median %.4f ms  (%.1f ns/eval)\n", med_c/1e6, med_c/nTotal)
    push!(perf_results, ("carpet_4x200", 800, med_c/nTotal,
          minimum(b_c).time/nTotal, maximum(b_c).time/nTotal,
          quantile(b_c.times, 0.05)/nTotal, quantile(b_c.times, 0.95)/nTotal))
    println()

    # ── B5. Combined carpet 4x441=1764 ───────────────────────────────────
    println("--- B5. Combined slip carpet (4 x 441 = 1764) ---")
    kv = vec([k for k in range(-0.3,0.3,length=21), _ in range(-0.3,0.3,length=21)])
    av = vec([a for _ in range(-0.3,0.3,length=21), a in range(-0.3,0.3,length=21)])
    Nc = 21*21; Nt = length(Fz_vals) * Nc
    mat_cb = Matrix{Float64}(undef, Nt, 6)
    row = 1
    for Fz in Fz_vals, j in 1:Nc
        mat_cb[row, :] = [Fz, kv[j], av[j], 0.0, 0.0, 16.7]
        row += 1
    end
    out_cb = Matrix{Float64}(undef, Nt, 30)
    b_cb = @benchmark mfeval!($out_cb, $p61, $mat_cb, $M111) evals=10 samples=100 seconds=5
    med_cb = median(b_cb).time
    @printf("  %d evals: median %.4f ms  (%.1f ns/eval)\n", Nt, med_cb/1e6, med_cb/Nt)
    push!(perf_results, ("carpet_combined", Nt, med_cb/Nt,
          minimum(b_cb).time/Nt, maximum(b_cb).time/Nt,
          quantile(b_cb.times, 0.05)/Nt, quantile(b_cb.times, 0.95)/Nt))
    println()

    return perf_results
end

perf_results = run_benchmarks()

# ── Write Julia performance CSV ───────────────────────────────────────────
println("--- Writing performance results to $PERF_CSV ---")
open(PERF_CSV, "w") do io
    println(io, "test_name,N,median_ns_per_eval,min_ns,max_ns,p05_ns,p95_ns")
    for (name, n, med, mn, mx, p05, p95) in perf_results
        @printf(io, "%s,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n", name, n, med, mn, mx, p05, p95)
    end
end

# ── Side-by-side comparison if MATLAB perf CSV exists ─────────────────────
let matlab_perf_path = joinpath(@__DIR__, "crossval_perf_matlab.csv")
    if isfile(matlab_perf_path)
        println("\n--- Side-by-side comparison ---\n")

        matlab_data = Dict{String, Float64}()
        for line in readlines(matlab_perf_path)[2:end]
            parts = split(line, ',')
            matlab_data[parts[1]] = parse(Float64, parts[3])  # median_us_per_eval
        end

        @printf("  %-25s  %12s  %12s  %10s\n", "Test", "MATLAB (us)", "Julia (ns)", "Speedup")
        @printf("  %-25s  %12s  %12s  %10s\n", "-"^25, "-"^12, "-"^12, "-"^10)

        open(COMP_CSV, "w") do io
            println(io, "test_name,matlab_us_per_eval,julia_ns_per_eval,speedup_x")
            for (name, _, med_ns, _, _, _, _) in perf_results
                if haskey(matlab_data, name)
                    matlab_us = matlab_data[name]
                    julia_us  = med_ns / 1000.0
                    speedup   = matlab_us / julia_us
                    @printf("  %-25s  %10.2f us  %10.1f ns  %8.1fx\n",
                            name, matlab_us, med_ns, speedup)
                    @printf(io, "%s,%.4f,%.2f,%.2f\n", name, matlab_us, med_ns, speedup)
                end
            end
        end
        println("\n  Comparison saved to $COMP_CSV")
    else
        println("\n  (No MATLAB perf CSV found — copy crossval_perf_matlab.csv to test/)")
    end
end

# ── Final summary ─────────────────────────────────────────────────────────
println("\n=============================================================")
@printf("  Validation: %d/%d tests passed", val_result.total_pass, val_result.total_tests)
if val_result.total_fail > 0
    @printf(" (%d FAILED)", val_result.total_fail)
end
println()
println("  Performance: $(length(perf_results)) benchmarks completed")
println("=============================================================")