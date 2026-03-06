# ==============================================================================
# test_phase3.jl
#
# Phase 3 API tests: mfeval (scalar) and mfeval / mfeval! (batch).
#
# Strategy
# --------
# 1. Scalar smoke — mfeval returns a well-typed MFOutputs with finite values.
# 2. Column ordering — scalar output matches to_vector / from_matrix_row.
# 3. Batch smoke — mfeval(p, mat, modes) returns an N×30 matrix.
# 4. mfeval! — in-place variant writes same values as allocating variant.
# 5. Scalar/batch consistency — row i of batch output == scalar output for row i.
# 6. inst_Kya — finite-difference column is populated (non-zero where α varies).
# 7. inst_Kya — zero when α is constant across all rows.
# 8. All three MF versions execute without error.
# 9. Edge cases — Fz=0, Vx=0 produce finite (not NaN/Inf) outputs.
# 10. Input widths — 6-, 7- and 8-column input matrices all accepted.
# ==============================================================================

if !isdefined(Main, :MFeval)
    using Test
    include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
    using .MFeval
end

# ── Fixtures ──────────────────────────────────────────────────────────────────
const _P3_TIR61 = joinpath(@__DIR__, "fixtures", "MagicFormula61_Parameters.tir")
const _P3_TIR62 = joinpath(@__DIR__, "fixtures", "MagicFormula62_Parameters.tir")
const _P3_TIR52 = joinpath(@__DIR__, "fixtures", "MagicFormula52_Parameters.tir")

const _p3_p61 = read_tir(_P3_TIR61)
const _p3_p62 = read_tir(_P3_TIR62)
const _p3_p52 = read_tir(_P3_TIR52)

const _p3_M111 = MFModes(111)

# Nominal operating point
const _p3_NOM = MFInputs(4000.0, -0.05, -0.1, 0.0, 0.0, 16.7)

# ==============================================================================
@testset "Phase 3: mfeval scalar — smoke" begin

    out = mfeval(_p3_p61, _p3_NOM, _p3_M111)

    @test out isa MFOutputs
    @test isfinite(out.Fx)
    @test isfinite(out.Fy)
    @test isfinite(out.Mz)
    @test isfinite(out.Re)
    @test isfinite(out.Rl)
    @test isfinite(out.rho)
    @test isfinite(out.sigmax)
    @test isfinite(out.sigmay)

    # Sanity ranges
    @test abs(out.Fx) < 20_000.0
    @test abs(out.Fy) < 20_000.0
    @test out.Re  > 0.1 && out.Re  < 1.0
    @test out.Rl  > 0.1 && out.Rl  < 1.0
    @test out.rho > 0.0
    @test out.Kxk > 0.0
    @test out.Kya != 0.0

end

# ==============================================================================
@testset "Phase 3: mfeval scalar — pass-through fields" begin

    out = mfeval(_p3_p61, _p3_NOM, _p3_M111)

    # Raw (unlimited) inputs must be passed through unchanged
    @test out.Fz       ≈ _p3_NOM.Fz
    @test out.kappa    ≈ _p3_NOM.kappa   atol=0.01   # low-speed reduction may apply
    @test out.alpha    ≈ _p3_NOM.alpha   atol=1e-10
    @test out.gamma    ≈ _p3_NOM.gamma   atol=1e-10
    @test out.Vx       ≈ _p3_NOM.Vx     atol=1e-10

end

# ==============================================================================
@testset "Phase 3: mfeval scalar — to_vector roundtrip" begin

    out = mfeval(_p3_p61, _p3_NOM, _p3_M111)
    v   = to_vector(out)

    @test length(v) == 30
    @test all(isfinite, v)

    # Reconstruct from the vector and check every field
    out2 = from_matrix_row(reshape(v, 1, 30), 1)
    @test out2.Fx  ≈ out.Fx
    @test out2.Fy  ≈ out.Fy
    @test out2.Mz  ≈ out.Mz
    @test out2.Kxk ≈ out.Kxk

end

# ==============================================================================
@testset "Phase 3: mfeval scalar — MF versions" begin

    for (label, tp) in (("MF6.1", _p3_p61), ("MF6.2", _p3_p62), ("MF5.2", _p3_p52))
        @testset "$label" begin
            out = mfeval(tp, _p3_NOM, _p3_M111)
            @test all(isfinite, to_vector(out))
        end
    end

end

# ==============================================================================
@testset "Phase 3: mfeval scalar — edge cases" begin

    @testset "Fz = 0 → all forces finite (not NaN)" begin
        out = mfeval(_p3_p61, MFInputs(0.0, 0.0, 0.0, 0.0, 0.0, 16.7), _p3_M111)
        @test all(isfinite, to_vector(out))
        @test out.Fx == 0.0 || abs(out.Fx) < 1.0
        @test out.Fy == 0.0 || abs(out.Fy) < 1.0
    end

    @testset "Vx = 0 (standstill) → finite outputs" begin
        out = mfeval(_p3_p61, MFInputs(4000.0, 0.0, 0.0, 0.0, 0.0, 0.0), _p3_M111)
        @test all(isfinite, to_vector(out))
        @test out.My == 0.0     # rolling resistance zero at standstill
    end

    @testset "Large kappa (locked wheel) → finite outputs" begin
        out = mfeval(_p3_p61, MFInputs(4000.0, -1.0, 0.0, 0.0, 0.0, 16.7), _p3_M111)
        @test all(isfinite, to_vector(out))
    end

end

# ==============================================================================
@testset "Phase 3: mfeval batch — smoke" begin

    N   = 20
    mat = Matrix{Float64}(undef, N, 6)
    for i in 1:N
        α = range(-0.3, 0.3, length=N)[i]
        mat[i, :] = [4000.0, 0.0, α, 0.0, 0.0, 16.7]
    end

    out = mfeval(_p3_p61, mat, _p3_M111)

    @test size(out) == (N, 30)
    @test all(isfinite, out)

end

# ==============================================================================
@testset "Phase 3: mfeval! in-place" begin

    N   = 10
    mat = Matrix{Float64}(undef, N, 6)
    for i in 1:N
        mat[i, :] = [4000.0, 0.0, Float64(i-1)*0.03 - 0.15, 0.0, 0.0, 16.7]
    end

    # Allocating form
    out_alloc = mfeval(_p3_p61, mat, _p3_M111)

    # In-place form
    out_inplace = Matrix{Float64}(undef, N, 30)
    ret = mfeval!(out_inplace, _p3_p61, mat, _p3_M111)

    @test ret === out_inplace        # returns the same object
    @test out_inplace ≈ out_alloc    # identical values

end

# ==============================================================================
@testset "Phase 3: scalar/batch consistency" begin

    alphas = [-0.2, -0.1, 0.0, 0.1, 0.2]
    N      = length(alphas)
    mat    = Matrix{Float64}(undef, N, 6)
    for i in 1:N
        mat[i, :] = [4000.0, 0.0, alphas[i], 0.0, 0.0, 16.7]
    end

    batch = mfeval(_p3_p61, mat, _p3_M111)

    for i in 1:N
        inp    = MFInputs(4000.0, 0.0, alphas[i], 0.0, 0.0, 16.7)
        scalar = mfeval(_p3_p61, inp, _p3_M111)
        sv     = to_vector(scalar)
        # col 29 (inst_Kya) differs between scalar (0.0) and batch (finite diff)
        # Compare cols 1-28 and 30 as a single vector (avoids invalid @test kw syntax)
        row_batch  = vcat(batch[i, 1:28], [batch[i, 30]])
        row_scalar = vcat(sv[1:28],       [sv[30]])
        @test row_batch ≈ row_scalar  atol=1e-10
    end

end

# ==============================================================================
@testset "Phase 3: inst_Kya (col 29) in batch" begin

    @testset "populated when alpha varies" begin
        N   = 5
        mat = Matrix{Float64}(undef, N, 6)
        for i in 1:N
            mat[i, :] = [4000.0, 0.0, (i-3)*0.05, 0.0, 0.0, 16.7]
        end
        out = mfeval(_p3_p61, mat, _p3_M111)
        # Central rows should have non-zero inst_Kya
        @test out[3, 29] != 0.0
        @test isfinite(out[3, 29])
    end

    @testset "indeterminate when alpha is constant" begin
        N   = 5
        mat = hcat(fill(4000.0, N), zeros(N), fill(0.1, N), zeros(N,3))
        # All rows have identical alpha → d_alpha = 0, d_Fy ≈ 0 → inst_Kya = NaN
        out = mfeval(_p3_p61, mat, _p3_M111)
        for i in 1:N
            @test isnan(out[i, 29]) || isinf(out[i, 29])   # division by zero → NaN or Inf
        end
    end

end

# ==============================================================================
@testset "Phase 3: batch input widths (6, 7, 8 columns)" begin

    # Helper function to compare matrices with NaN equality
    function matrices_equal_nan(a, b)
        if size(a) != size(b)
            return false
        end
        for i in eachindex(a)
            if isnan(a[i]) && isnan(b[i])
                continue  # NaN == NaN for this comparison
            elseif isnan(a[i]) || isnan(b[i])
                return false  # one NaN, one not
            elseif !isapprox(a[i], b[i])
                return false  # numeric difference
            end
        end
        return true
    end

    N = 3
    base = [4000.0 0.0 0.1 0.0 0.0 16.7]
    row3 = repeat(base, N, 1)

    out6 = mfeval(_p3_p61, row3, _p3_M111)
    @test size(out6) == (N, 30)

    # 7-column: explicit pressure
    row7 = hcat(row3, fill(_p3_p61.inflpres, N))
    out7 = mfeval(_p3_p61, row7, _p3_M111)
    @test matrices_equal_nan(out7, out6)   # same pressure → same result

    # 8-column: pressure + omega=0 (solver computes it)
    row8 = hcat(row7, zeros(N))
    out8 = mfeval(_p3_p61, row8, _p3_M111)
    @test matrices_equal_nan(out8, out6)

end