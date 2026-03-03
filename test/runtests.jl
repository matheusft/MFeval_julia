# ==============================================================================
# runtests.jl — MFeval.jl test suite entry point
#
# Run all tests with:
#   julia --project=. test/runtests.jl        (from the MFeval.jl/ root)
#   julia -e 'using Pkg; Pkg.test()'          (from anywhere, after Pkg.activate)
#
# Individual test files can also be run directly:
#   julia --project=. test/test_phase1.jl
# ==============================================================================

using Test

# ── Phase 1: types, I/O, structs ──────────────────────────────────────────────
println("=" ^ 60)
println("Phase 1 — Types, I/O and structs")
println("=" ^ 60)
include("test_phase1.jl")

# ── Future phases — uncomment as they are implemented ─────────────────────────
# println("=" ^ 60)
# println("Phase 2 — Solver kernels")
# println("=" ^ 60)
# include("test_mf61.jl")
# include("test_mf52.jl")
# include("test_mf62.jl")
# include("test_low_speed.jl")
# include("test_turn_slip.jl")

# println("=" ^ 60)
# println("Phase 4 — Numerical regression vs MATLAB")
# println("=" ^ 60)
# include("test_vs_matlab.jl")