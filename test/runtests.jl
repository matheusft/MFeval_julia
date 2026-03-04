# ==============================================================================
# runtests.jl — MFeval.jl test suite entry point
#
# Run all tests with:
#   julia --project=. test/runtests.jl        (from the MFeval.jl/ root)
#   julia -e 'using Pkg; Pkg.test()'          (from anywhere, after Pkg.activate)
#
# Individual test files can also be run directly (they self-bootstrap):
#   julia --project=. test/test_phase1.jl
#   julia --project=. test/test_phase2.jl
# ==============================================================================

using Test

# Load the module exactly once.  Both test files detect this and skip their
# own bootstrap when already included via runtests.jl.
include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
using .MFeval

# ── Phase 1: types, I/O, structs ──────────────────────────────────────────────
println("=" ^ 60)
println("Phase 1 — Types, I/O and structs")
println("=" ^ 60)
include(joinpath(@__DIR__, "test_phase1.jl"))

# ── Phase 2: scalar solver kernels ────────────────────────────────────────────
println("=" ^ 60)
println("Phase 2 — Scalar solver kernels")
println("=" ^ 60)
include(joinpath(@__DIR__, "test_phase2.jl"))

# ── Phase 3: public API ────────────────────────────────────────────────────────
println("=" ^ 60)
println("Phase 3 — Public API")
println("=" ^ 60)
include(joinpath(@__DIR__, "test_phase3.jl"))

# ── Phase 4: validation — edge cases, invariants, regression ──────────────────
println("=" ^ 60)
println("Phase 4 — Validation")
println("=" ^ 60)
include(joinpath(@__DIR__, "test_phase4.jl"))

nothing  # suppress REPL echo of the last @testset return value