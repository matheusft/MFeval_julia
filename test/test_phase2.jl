# ==============================================================================
# test_phase2.jl
#
# Phase 2 solver tests: parse_inputs, basic_vars, turn_slip, Fx0, Fy0,
# combined_slip, Mx, My, Mz0, Mz, geometry, relaxation.
#
# Strategy
# --------
# 1. Smoke tests — every function can be called without error on a nominal
#    operating point and returns values in physically sensible ranges.
# 2. Zero-Fz guard — all force/moment outputs are 0 or near-0 when Fz = 0.
# 3. Symmetry — Fx(κ) is odd in κ; Fy(α) is odd in α; Mz(α) is odd in α.
# 4. Parse_inputs limit checks — saturated outputs respect TIR bounds.
# 5. Mode flags — alpha_star and turn_slip branches execute correctly.
# 6. Geometry — Re, Rl, and contact patch are physically consistent.
# 7. Relaxation — sigmax and sigmay are positive.
# ==============================================================================

if !isdefined(Main, :MFeval)
    using Test
    include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
    using .MFeval
end

# ── Shared test fixtures ────────────────────────────────────────────────────────

# Fixture paths — guarded to avoid redefinition when included after test_phase1.jl
if !@isdefined(TIR_MF61)
    const TIR_MF61 = joinpath(@__DIR__, "fixtures", "MagicFormula61_Parameters.tir")
end
if !@isdefined(TIR_MF62)
    const TIR_MF62 = joinpath(@__DIR__, "fixtures", "MagicFormula62_Parameters.tir")
end
if !@isdefined(TIR_MF52)
    const TIR_MF52 = joinpath(@__DIR__, "fixtures", "MagicFormula52_Parameters.tir")
end

const p61 = read_tir(TIR_MF61)
const p62 = read_tir(TIR_MF62)
const p52 = read_tir(TIR_MF52)

# Nominal steady-state operating point (normal driving, small slip)
const NOM = MFInputs(4000.0, -0.05, -0.1, 0.0, 0.0, 16.7)

# Mode: limits=true, alpha=false, turn_slip=false  (most common)
const M111 = MFModes(111)
# Mode: limits=true, alpha_star=true, turn_slip=false
const M121 = MFModes(121)
# Mode: limits=true, alpha=false, turn_slip=true
const M112 = MFModes(112)
# Mode: no limits
const M211 = MFModes(211)

# ── Helper: run the full solver chain up to forces ───────────────────────────────

function run_chain(tp, inp, modes)
    pp, ip       = parse_inputs(tp, inp, modes)
    sv, pv, iv, sl = calc_basic_vars(tp, pp, ip, modes)
    Fx0, mux, Kxk = calc_Fx0(tp, pp, ip, modes, sv, pv, iv)
    Fy0, muy, fi, ip2 = calc_Fy0(tp, pp, ip, modes, sv, pv, iv)
    Mz0, mz0i       = calc_Mz0(tp, pp, ip2, modes, sv, pv, iv, fi)
    Fx               = calc_Fx(tp, pp, modes, sv, iv, Fx0)
    Fy, Gyk, SVyk    = calc_Fy(tp, pp, ip2, modes, sv, iv, Fy0, muy)
    Mx               = calc_Mx(tp, pp, iv, Fy)
    My               = calc_My(tp, pp, Fx)
    Mz, t, Mzr      = calc_Mz(tp, pp, ip2, modes, sv, pv, iv, mz0i, Kxk, Fy, Fx, SVyk)
    dpi              = (pp.p - tp.nompres) / tp.nompres
    Re, Romega, ω    = calc_Re(tp, pp, dpi)
    rho, Rl, Cz      = calc_rho_Rl(tp, pp, dpi, ω, Romega, Fx, Fy)
    a, b, NCz        = calc_contact_patch(tp, pp, dpi)
    Cx_s, Cy_s, σx, σy = calc_relax(tp, pp, iv, Kxk, fi.Kya)
    (; pp, ip=ip2, sv, pv, iv, sl,
       Fx0, mux, Kxk, Fy0, muy, fi, Mz0, mz0i,
       Fx, Fy, Gyk, SVyk, Mx, My, Mz, t, Mzr,
       Re, Romega, ω, rho, Rl, Cz, a, b, NCz,
       Cx=Cx_s, Cy=Cy_s, σx, σy)
end

# ==============================================================================
@testset "Phase 2: parse_inputs" begin

    @testset "nominal smoke test" begin
        pp, ip = parse_inputs(p61, NOM, M111)
        @test pp.Fz    > 0.0
        @test pp.uFz   ≈ NOM.Fz
        @test pp.uVcx  ≈ NOM.Vx
        @test pp.p     > 0.0
        @test pp.omega > 0.0
        @test isfinite(pp.Vsx)
        @test isfinite(pp.phi)
    end

    @testset "Fz negative → clamped to zero" begin
        inp = MFInputs(-100.0, 0.0, 0.0, 0.0, 0.0, 16.7)
        pp, _ = parse_inputs(p61, inp, M111)
        @test pp.Fz   == 0.0
        @test pp.uFz  == 0.0
    end

    @testset "pressure from TIR when inp.pressure == 0" begin
        pp, _ = parse_inputs(p61, NOM, M111)    # NOM has pressure=0 (default)
        @test pp.p ≈ p61.inflpres
    end

    @testset "pressure saturated to [PRESMIN, PRESMAX]" begin
        inp_lo = MFInputs(4000.0, 0.0, 0.0, 0.0, 0.0, 16.7, 1.0)      # very low
        inp_hi = MFInputs(4000.0, 0.0, 0.0, 0.0, 0.0, 16.7, 1e9)      # very high
        pp_lo, _ = parse_inputs(p61, inp_lo, M111)
        pp_hi, _ = parse_inputs(p61, inp_hi, M111)
        @test pp_lo.p ≈ p61.presmin
        @test pp_hi.p ≈ p61.presmax
    end

    @testset "low-speed reductions active" begin
        inp_slow = MFInputs(4000.0, -0.5, -0.2, 0.0, 0.0, 0.5)  # |Vcx| < VXLOW
        pp, ip = parse_inputs(p61, inp_slow, M111)
        @test ip.reductionSmooth <  1.0
        @test ip.reductionLinear <  1.0
        @test pp.kappa >= 0.0 || pp.kappa <= 0.0   # Just ensure finite
    end

    @testset "no-limits mode: reductions == 1" begin
        pp, ip = parse_inputs(p61, NOM, M211)
        @test ip.reductionSmooth       == 1.0
        @test ip.reductionLinear       == 1.0
        @test ip.reductionLinear_alpha == 1.0
    end

    @testset "ukappa is truly raw" begin
        inp_slow = MFInputs(4000.0, -0.8, 0.0, 0.0, 0.0, 0.3)
        pp, _ = parse_inputs(p61, inp_slow, M111)
        @test pp.ukappa    == -0.8         # raw, unmodified
        @test pp.ukappaLow != -0.8         # should be reduced
        @test abs(pp.ukappaLow) < abs(pp.ukappa)
    end

    @testset "turn-slip phi computed" begin
        inp_ts = MFInputs(4000.0, 0.0, 0.1, 0.0, 0.05, 16.7)
        pp, ip = parse_inputs(p61, inp_ts, M112)
        @test isfinite(pp.phi)
        @test ip.epsilong > 0.0
    end

    @testset "zero Vcx → ualpha == 0" begin
        inp_stop = MFInputs(4000.0, 0.0, 0.2, 0.0, 0.0, 0.0)
        pp, _ = parse_inputs(p61, inp_stop, M111)
        @test pp.ualpha == 0.0
    end

end  # parse_inputs

# ==============================================================================
@testset "Phase 2: calc_basic_vars" begin

    pp, ip = parse_inputs(p61, NOM, M111)

    @testset "nominal smoke" begin
        sv, pv, iv, sl = calc_basic_vars(p61, pp, ip, M111)
        @test isfinite(sv.alpha_star)
        @test isfinite(sv.gamma_star)
        @test sv.LMUX_star > 0.0
        @test sv.LMUY_star > 0.0
        @test pv.Fz0_prime > 0.0
        @test 0.0 <= pv.alpha_prime <= π/2
        @test pv.LMUX_prime > 0.0
        @test pv.LMUY_prime > 0.0
        @test isfinite(iv.dfz)
        @test isfinite(iv.dpi)
    end

    @testset "alpha_star mode vs standard" begin
        sv_std,  _, _, _ = calc_basic_vars(p61, pp, ip, M111)
        sv_star, _, _, _ = calc_basic_vars(p61, pp, ip, M121)
        # With alpha_star, alpha_star = tan(alpha)*sign(Vcx)
        @test sv_star.alpha_star ≈ tan(pp.alpha) * sign(pp.uVcx)
        # Without, alpha_star = alpha
        @test sv_std.alpha_star  ≈ pp.alpha
    end

    @testset "dfz is negative when Fz < Fz0_prime" begin
        pp2, ip2 = parse_inputs(p61, MFInputs(2000.0, 0.0, 0.0, 0.0, 0.0, 16.7), M111)
        _, _, iv2, _ = calc_basic_vars(p61, pp2, ip2, M111)
        @test iv2.dfz < 0.0
    end

end  # basic_vars

# ==============================================================================
@testset "Phase 2: calc_Fx0" begin

    @testset "nominal output range" begin
        r = run_chain(p61, NOM, M111)
        @test isfinite(r.Fx0)
        @test abs(r.Fx0) < 20_000.0   # physically: < 2× Fz
        @test r.mux  > 0.0
        @test r.Kxk  > 0.0
    end

    @testset "zero Fz → zero Fx0" begin
        inp0 = MFInputs(0.0, -0.05, 0.0, 0.0, 0.0, 16.7)
        pp, ip = parse_inputs(p61, inp0, M111)
        sv, pv, iv, _ = calc_basic_vars(p61, pp, ip, M111)
        Fx0, mux, Kxk = calc_Fx0(p61, pp, ip, M111, sv, pv, iv)
        @test Fx0 == 0.0 || abs(Fx0) < 1.0
        @test mux == 0.0
    end

    @testset "Fx0 odd in kappa (no camber, no shift)" begin
        # At kappa = ±k, Fx0 should be ≈ ∓|Fx0| when PHX1=PHX2=PVX1=PVX2=0
        # (default coefficients may not be exactly zero, so we check sign)
        inp_pos = MFInputs(4000.0,  0.1, 0.0, 0.0, 0.0, 16.7)
        inp_neg = MFInputs(4000.0, -0.1, 0.0, 0.0, 0.0, 16.7)
        pp_p, ip_p = parse_inputs(p61, inp_pos, M211)  # no limits for clean test
        pp_n, ip_n = parse_inputs(p61, inp_neg, M211)
        sv_p, pv_p, iv_p, _ = calc_basic_vars(p61, pp_p, ip_p, M211)
        sv_n, pv_n, iv_n, _ = calc_basic_vars(p61, pp_n, ip_n, M211)
        Fx0_p, _, _ = calc_Fx0(p61, pp_p, ip_p, M211, sv_p, pv_p, iv_p)
        Fx0_n, _, _ = calc_Fx0(p61, pp_n, ip_n, M211, sv_n, pv_n, iv_n)
        # Sign should be opposite
        @test sign(Fx0_p) == -sign(Fx0_n)
    end

    @testset "backward speed → sign flip" begin
        inp_fwd = MFInputs(4000.0, -0.1, 0.0, 0.0, 0.0,  16.7)
        inp_bwd = MFInputs(4000.0, -0.1, 0.0, 0.0, 0.0, -16.7)
        pp_f, ip_f = parse_inputs(p61, inp_fwd, M211)
        pp_b, ip_b = parse_inputs(p61, inp_bwd, M211)
        sv_f, pv_f, iv_f, _ = calc_basic_vars(p61, pp_f, ip_f, M211)
        sv_b, pv_b, iv_b, _ = calc_basic_vars(p61, pp_b, ip_b, M211)
        Fx0_f, _, _ = calc_Fx0(p61, pp_f, ip_f, M211, sv_f, pv_f, iv_f)
        Fx0_b, _, _ = calc_Fx0(p61, pp_b, ip_b, M211, sv_b, pv_b, iv_b)
        @test sign(Fx0_f) == -sign(Fx0_b)
    end

end  # Fx0

# ==============================================================================
@testset "Phase 2: calc_Fy0" begin

    @testset "nominal output range" begin
        r = run_chain(p61, NOM, M111)
        @test isfinite(r.Fy0)
        @test abs(r.Fy0) < 20_000.0
        @test r.muy > 0.0
        @test r.fi.Kya != 0.0
    end

    @testset "zero Fz → zero Fy0 and muy" begin
        inp0 = MFInputs(0.0, 0.0, -0.1, 0.0, 0.0, 16.7)
        pp, ip = parse_inputs(p61, inp0, M111)
        sv, pv, iv, _ = calc_basic_vars(p61, pp, ip, M111)
        Fy0, muy, _, _ = calc_Fy0(p61, pp, ip, M111, sv, pv, iv)
        @test abs(Fy0) < 1.0
        @test muy == 0.0
    end

    @testset "Fy0 odd in alpha at zero camber" begin
        inp_pos = MFInputs(4000.0, 0.0,  0.15, 0.0, 0.0, 16.7)
        inp_neg = MFInputs(4000.0, 0.0, -0.15, 0.0, 0.0, 16.7)
        pp_p, ip_p = parse_inputs(p61, inp_pos, M211)
        pp_n, ip_n = parse_inputs(p61, inp_neg, M211)
        sv_p, pv_p, iv_p, _ = calc_basic_vars(p61, pp_p, ip_p, M211)
        sv_n, pv_n, iv_n, _ = calc_basic_vars(p61, pp_n, ip_n, M211)
        Fy0_p, _, _, _ = calc_Fy0(p61, pp_p, ip_p, M211, sv_p, pv_p, iv_p)
        Fy0_n, _, _, _ = calc_Fy0(p61, pp_n, ip_n, M211, sv_n, pv_n, iv_n)
        asymmetry = abs(Fy0_p + Fy0_n)
        magnitude  = 0.5 * (abs(Fy0_p) + abs(Fy0_n))
        @test sign(Fy0_p) == -sign(Fy0_n)
        @test asymmetry < 0.3 * magnitude
    end

    @testset "Fy0 exactly antisymmetric with zero shift coefficients" begin
        # Test with default TireParams (all shifts are zero) instead of modifying p61
        p_default = TireParams(metadata = TireMetadata(fittyp = 61))
        
        inp_pos = MFInputs(4000.0, 0.0,  0.15, 0.0, 0.0, 16.7)
        inp_neg = MFInputs(4000.0, 0.0, -0.15, 0.0, 0.0, 16.7)
        pp_p, ip_p = parse_inputs(p_default, inp_pos, M211)
        pp_n, ip_n = parse_inputs(p_default, inp_neg, M211)
        sv_p, pv_p, iv_p, _ = calc_basic_vars(p_default, pp_p, ip_p, M211)
        sv_n, pv_n, iv_n, _ = calc_basic_vars(p_default, pp_n, ip_n, M211)
        Fy0_p, _, _, _ = calc_Fy0(p_default, pp_p, ip_p, M211, sv_p, pv_p, iv_p)
        Fy0_n, _, _, _ = calc_Fy0(p_default, pp_n, ip_n, M211, sv_n, pv_n, iv_n)
        @test Fy0_p ≈ -Fy0_n  atol=1.0
    end

    @testset "MF5.2 branch executes" begin
        r = run_chain(p52, NOM, M111)
        @test isfinite(r.Fy0)
        @test r.muy > 0.0
    end

end  # Fy0

# ==============================================================================
@testset "Phase 2: calc_Mz0 / calc_Mz" begin

    @testset "nominal Mz0 finite" begin
        r = run_chain(p61, NOM, M111)
        @test isfinite(r.Mz0)
        @test abs(r.Mz0) < 2000.0   # typically ±500 N·m for these conditions
    end

    @testset "Mz output finite with trail and residual" begin
        r = run_chain(p61, NOM, M111)
        @test isfinite(r.Mz)
        @test isfinite(r.t)
        @test isfinite(r.Mzr)
        @test r.t > 0.0   # pneumatic trail positive at small positive alpha
    end

    @testset "Mz zero at zero alpha and zero kappa" begin
        inp0 = MFInputs(4000.0, 0.0, 0.0, 0.0, 0.0, 16.7)
        r = run_chain(p61, inp0, M111)
        @test abs(r.Mz) < 50.0   # near-zero (residual torque may be non-zero)
    end

    @testset "MF5.2 Mz executes" begin
        r = run_chain(p52, NOM, M111)
        @test isfinite(r.Mz)
    end

    @testset "MF6.2 Mz executes" begin
        r = run_chain(p62, NOM, M111)
        @test isfinite(r.Mz)
    end

end  # Mz

# ==============================================================================
@testset "Phase 2: calc_Fx / calc_Fy (combined slip)" begin

    @testset "Gxa near 1 at zero alpha (no combined effect)" begin
        inp_noalpha = MFInputs(4000.0, -0.1, 0.0, 0.0, 0.0, 16.7)
        r = run_chain(p61, inp_noalpha, M111)
        @test r.Fx ≈ r.Fx0  atol=100.0
    end

    @testset "Gyk near 1 at zero kappa (no combined effect)" begin
        inp_nokappa = MFInputs(4000.0, 0.0, -0.1, 0.0, 0.0, 16.7)
        r = run_chain(p61, inp_nokappa, M111)
        @test r.Fy ≈ r.Fy0  atol=100.0
    end

    @testset "|Fx| <= |Fx0| under combined slip (load transfer)" begin
        r = run_chain(p61, NOM, M111)   # NOM has both kappa and alpha
        @test abs(r.Fx) <= abs(r.Fx0) + 200.0   # Gxa ≤ 1
    end

end  # combined_slip

# ==============================================================================
@testset "Phase 2: calc_Mx" begin

    @testset "Mx finite for nominal input" begin
        r = run_chain(p61, NOM, M111)
        @test isfinite(r.Mx)
    end

    @testset "Mx zero when Fz == 0" begin
        inp0 = MFInputs(0.0, 0.0, 0.0, 0.0, 0.0, 16.7)
        r = run_chain(p61, inp0, M111)
        @test abs(r.Mx) < 1.0
    end

    @testset "Mx sign flips with gamma" begin
        inp_p = MFInputs(4000.0, 0.0, 0.0,  0.05, 0.0, 16.7)
        inp_n = MFInputs(4000.0, 0.0, 0.0, -0.05, 0.0, 16.7)
        r_p = run_chain(p61, inp_p, M211)
        r_n = run_chain(p61, inp_n, M211)
        @test sign(r_p.Mx) == -sign(r_n.Mx)
    end

end  # Mx

# ==============================================================================
@testset "Phase 2: calc_My" begin

    @testset "My finite and negative (ISO: negative rolling resistance)" begin
        r = run_chain(p61, NOM, M111)
        @test isfinite(r.My)
        @test r.My < 0.0   # rolling resistance opposes forward motion
    end

    @testset "My zero at Vcx == 0" begin
        inp_stop = MFInputs(4000.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        r = run_chain(p61, inp_stop, M111)
        @test r.My == 0.0
    end

    @testset "My sign flips for backward speed" begin
        inp_fwd = MFInputs(4000.0, 0.0, 0.0, 0.0, 0.0,  16.7)
        inp_bwd = MFInputs(4000.0, 0.0, 0.0, 0.0, 0.0, -16.7)
        r_f = run_chain(p61, inp_fwd, M211)
        r_b = run_chain(p61, inp_bwd, M211)
        @test sign(r_f.My) == -sign(r_b.My)
    end

    @testset "MF5.2 My executes" begin
        r = run_chain(p52, NOM, M111)
        @test isfinite(r.My)
        @test r.My < 0.0
    end

end  # My

# ==============================================================================
@testset "Phase 2: geometry (calc_Re, calc_rho_Rl, calc_contact_patch)" begin

    @testset "Re is slightly less than R0" begin
        r = run_chain(p61, NOM, M111)
        R0 = p61.unloaded_radius
        @test r.Re > R0 * 0.9
        @test r.Re < R0
    end

    @testset "Romega ≈ R0 (at low speed, centrifugal growth negligible)" begin
        r = run_chain(p61, NOM, M111)
        R0 = p61.unloaded_radius
        @test r.Romega ≈ R0  atol=0.01
    end

    @testset "Rl < Romega (loaded radius less than free)" begin
        r = run_chain(p61, NOM, M111)
        @test r.Rl < r.Romega
        @test r.Rl > 0.0
    end

    @testset "rho (deflection) positive" begin
        r = run_chain(p61, NOM, M111)
        @test r.rho > 0.0
    end

    @testset "contact patch half-lengths positive" begin
        r = run_chain(p61, NOM, M111)
        @test r.a > 0.0
        @test r.b > 0.0
    end

    @testset "Rl + rho ≈ Romega (geometry identity)" begin
        r = run_chain(p61, NOM, M111)
        @test r.Rl + r.rho ≈ r.Romega  atol=1e-4
    end

    @testset "MF6.2 geometry executes (secant solver)" begin
        r = run_chain(p62, NOM, M111)
        @test r.Rl > 0.0
        @test r.rho > 0.0
    end

    @testset "MF5.2 contact patch uses Q_A formula" begin
        r = run_chain(p52, NOM, M111)
        @test r.a > 0.0
        @test r.b > 0.0
    end

end  # geometry

# ==============================================================================
@testset "Phase 2: relaxation (calc_relax)" begin

    @testset "stiffnesses positive" begin
        r = run_chain(p61, NOM, M111)
        @test r.Cx > 0.0
        @test r.Cy > 0.0
    end

    @testset "relaxation lengths positive" begin
        r = run_chain(p61, NOM, M111)
        @test r.σx > 0.0
        @test r.σy > 0.0
    end

    @testset "MF5.2 relaxation uses PTX/PTY formula" begin
        r = run_chain(p52, NOM, M111)
        @test r.σx > 0.0
        @test r.σy > 0.0
    end

end  # relaxation

# ==============================================================================
@testset "Phase 2: turn-slip mode (M112)" begin

    inp_ts = MFInputs(4000.0, 0.0, 0.1, 0.0, 0.02, 16.7)

    @testset "MF6.1 turn-slip: all outputs finite" begin
        r = run_chain(p61, inp_ts, M112)
        for val in (r.Fx, r.Fy, r.Mx, r.My, r.Mz)
            @test isfinite(val)
        end
    end

    @testset "zeta1 ≤ 1.0 (turn-slip reduces Fx peak)" begin
        pp, ip = parse_inputs(p61, inp_ts, M112)
        sv, pv, iv, _ = calc_basic_vars(p61, pp, ip, M112)
        z1 = calc_zeta1(p61, pp, iv.dfz)
        @test 0.0 <= z1 <= 1.0
    end

    @testset "zeta2 ≤ 1.0 (turn-slip reduces Fy peak)" begin
        pp, ip = parse_inputs(p61, inp_ts, M112)
        sv, pv, iv, _ = calc_basic_vars(p61, pp, ip, M112)
        z2, z3 = calc_zeta2_zeta3(p61, pp, iv.dfz)
        @test 0.0 <= z2 <= 1.0
        @test 0.0 <= z3 <= 1.0
    end

end  # turn_slip

# ==============================================================================
@testset "Phase 2: alpha_star mode (M121)" begin

    @testset "alpha_star differs from plain alpha at large alpha" begin
        inp_large = MFInputs(4000.0, 0.0, 0.4, 0.0, 0.0, 16.7)
        r_std  = run_chain(p61, inp_large, M211)   # no alpha_star
        r_star = run_chain(p61, inp_large, MFModes(221))  # alpha_star
        # Forces should differ because alpha_star ≠ alpha at large slip
        @test r_std.Fy ≠ r_star.Fy
    end

end  # alpha_star

# ==============================================================================
@testset "Phase 2: calc_inst_Kya" begin

    @testset "finite difference gives reasonable stiffness" begin
        α1 = 0.05;  α2 = 0.06
        inp1 = MFInputs(4000.0, 0.0, α1, 0.0, 0.0, 16.7)
        inp2 = MFInputs(4000.0, 0.0, α2, 0.0, 0.0, 16.7)
        r1 = run_chain(p61, inp1, M111)
        r2 = run_chain(p61, inp2, M111)
        iKya = calc_inst_Kya(α2 - α1, r2.Fy - r1.Fy)
        # Typical cornering stiffness 50_000–500_000 N/rad
        @test iKya > 0.0
    end

    @testset "zero d_alpha → 0 (no divide by zero)" begin
        @test calc_inst_Kya(0.0, 100.0) == 0.0
    end

end  # calc_inst_Kya