# ==============================================================================
# test_phase1.jl
#
# Phase 1 tests: MFModes, MFInputs, MFOutputs, TireParams, read_tir.
#
# When run via runtests.jl the module is already loaded.
# When run standalone the bootstrap block below loads it directly.
# ==============================================================================

if !isdefined(Main, :MFeval)
    using Test
    include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
    using .MFeval
end

# Fixture TIR files live next to this test file, under test/fixtures/
const TIR_DIR = joinpath(@__DIR__, "fixtures")
const TIR_MF52 = joinpath(TIR_DIR, "MagicFormula52_Parameters.tir")
const TIR_MF61 = joinpath(TIR_DIR, "MagicFormula61_Parameters.tir")
const TIR_MF62 = joinpath(TIR_DIR, "MagicFormula62_Parameters.tir")

# ==============================================================================
@testset "MFModes decoding" begin
    m = MFModes(111)
    @test m.use_limits_check == true
    @test m.use_alpha_star   == false
    @test m.use_turn_slip    == false
    @test m.is_valid         == true

    m = MFModes(121)
    @test m.use_limits_check == true
    @test m.use_alpha_star   == true
    @test m.use_turn_slip    == false
    @test m.is_valid         == true

    m = MFModes(112)
    @test m.use_limits_check == true
    @test m.use_alpha_star   == false
    @test m.use_turn_slip    == true
    @test m.is_valid         == true

    m = MFModes(222)
    @test m.use_limits_check == false
    @test m.use_alpha_star   == true
    @test m.use_turn_slip    == true
    @test m.is_valid         == true

    # Invalid useModes are now accepted (matches MATLAB) but marked as invalid
    m_invalid1 = MFModes(311)  # bad hundreds digit
    @test m_invalid1.is_valid == false
    
    m_invalid2 = MFModes(131)  # bad tens digit
    @test m_invalid2.is_valid == false
    
    m_invalid3 = MFModes(113)  # bad units digit
    @test m_invalid3.is_valid == false
end

# ==============================================================================
@testset "MFInputs constructors" begin
    inp6 = MFInputs(4000.0, 0.0, 0.05, 0.0, 0.0, 16.0)
    @test inp6.Fz       == 4000.0
    @test inp6.kappa    == 0.0
    @test inp6.alpha    == 0.05
    @test inp6.gamma    == 0.0
    @test inp6.phit     == 0.0
    @test inp6.Vx       == 16.0
    @test inp6.pressure == 0.0
    @test inp6.omega    == 0.0

    inp7 = MFInputs(4000.0, 0.0, 0.05, 0.0, 0.0, 16.0, 250000.0)
    @test inp7.pressure == 250000.0
    @test inp7.omega    == 0.0

    inp8 = MFInputs(4000.0, 0.0, 0.05, 0.0, 0.0, 16.0, 250000.0, 45.0)
    @test inp8.pressure == 250000.0
    @test inp8.omega    == 45.0
end

# ==============================================================================
@testset "MFOutputs round-trip" begin
    v = collect(Float64, 1:30)
    m = reshape(v, 1, 30)
    o = from_matrix_row(m, 1)
    @test to_vector(o) == v
end

# ==============================================================================
@testset "TireParams default construction" begin
    meta = TireMetadata(fittyp=61)
    p    = TireParams(metadata=meta)

    @test is_mf61(p)
    @test version(p) === :MF61

    @test p.longvl           ≈ 16.7
    @test p.unloaded_radius  ≈ 0.3135
    @test p.fnomin           ≈ 4000.0
    @test p.nompres          ≈ 220000.0
    @test p.lfzo             ≈ 1.0
    @test p.lcx              ≈ 1.0
    @test p.lmux             ≈ 1.0
    @test p.pcx1             ≈ 1.65
    @test p.pdx1             ≈ 1.3
    @test p.pky1             ≈ -20.0
    @test p.qsy7             ≈ 0.85
    @test p.qsy8             ≈ -0.4
    @test p.presmin          ≈ 10_000.0
    @test p.presmax          ≈ 1_000_000.0

    meta52 = TireMetadata(fittyp=21)
    p52 = TireParams(metadata=meta52)
    @test is_mf52(p52)
    @test version(p52) === :MF52

    meta62 = TireMetadata(fittyp=62)
    p62 = TireParams(metadata=meta62)
    @test is_mf62(p62)
    @test version(p62) === :MF62

    meta_bad = TireMetadata(fittyp=99)
    @test_throws Exception TireParams(metadata=meta_bad)
end

# ==============================================================================
@testset "read_tir MF6.1" begin
    if !isfile(TIR_MF61)
        @warn "MF6.1 TIR file not found at $TIR_MF61, skipping file-based tests"
    else
        p = read_tir(TIR_MF61)

        @test is_mf61(p)
        @test p.metadata.fittyp == 61

        @test p.unloaded_radius  ≈ 0.3943      atol=1e-6
        @test p.fnomin           ≈ 6752.0      atol=1e-3
        @test p.nompres          ≈ 260000.0    atol=1e-1

        @test p.lfzo  ≈ 1.0
        @test p.lcx   ≈ 1.0
        @test p.lmux  ≈ 1.0
        @test p.lmuy  ≈ 1.0

        @test p.pcx1  ≈ 1.5482   atol=1e-4
        @test p.pdx1  ≈ 1.1632   atol=1e-4
        @test p.pkx1  ≈ 32.9102  atol=1e-4

        @test p.ppx1  ≈ -0.8733  atol=1e-4
        @test p.ppy1  ≈  0.34039 atol=1e-4

        @test p.pdxp1 ≈ 0.4      atol=1e-4
        @test p.pecp1 ≈ 0.5      atol=1e-4

        @test p.kpumin ≈ -1.0
        @test p.kpumax ≈  1.0
    end
end

# ==============================================================================
@testset "read_tir MF5.2" begin
    if !isfile(TIR_MF52)
        @warn "MF5.2 TIR file not found at $TIR_MF52, skipping"
    else
        p = read_tir(TIR_MF52)
        @test is_mf52(p)
        @test p.metadata.fittyp ∈ (6, 21)
    end
end

# ==============================================================================
@testset "read_tir MF6.2" begin
    if !isfile(TIR_MF62)
        @warn "MF6.2 TIR file not found at $TIR_MF62, skipping"
    else
        p = read_tir(TIR_MF62)
        @test is_mf62(p)
        @test p.metadata.fittyp == 62
        @test p.q_cam1 isa Float64
    end
end

# ==============================================================================
@testset "read_tir error handling" begin
    @test_throws Exception read_tir("nonexistent_file.tir")
end

# ==============================================================================
@testset "TireParams type concreteness" begin
    meta = TireMetadata(fittyp=61)
    p    = TireParams(metadata=meta)
    for fname in fieldnames(TireParams)
        fname == :metadata && continue
        @test fieldtype(TireParams{:MF61}, fname) === Float64
    end
end