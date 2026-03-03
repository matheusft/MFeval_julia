# ==============================================================================
# My.jl
#
# Rolling resistance moment My.
# Mirrors Solver.calculateMy from the reference MATLAB code.
#
# The book equation (4.E70) is NOT used because it gives wrong sign and
# has an error (Fz instead of Fz0 in the first term). Paper Eqn (A48)
# is used for MF6.1/6.2; the MF-Tyre manual equation for MF5.2.
#
# Reference: Pacejka (2012) Book Appendix 3, Eqn (A48)
#            MF-Tyre equation manual (MF5.2)
# ==============================================================================

"""
    calc_My(p, pp, Fx) → My

Rolling resistance moment (N·m).

Applies:
- Low-Fz empirical correction  (Fz_eff = Fz·(Fz/FZMIN) when Fz < FZMIN)
- Low-speed interpolation to smooth out the zero-crossing near standstill
- Backward-speed sign flip
- Zero-speed zeroing
"""
@inline @fastmath function calc_My(p  ::TireParams,
                                    pp ::PostProInputs,
                                    Fx ::Float64) ::Float64

    Vcx   = pp.uVcx
    kappa = pp.ukappa
    gamma = pp.gamma
    press = pp.p
    uFz   = pp.uFz

    # Low-Fz empirical correction (matches reference MATLAB)
    Fz_eff = (uFz < p.fzmin) ? uFz * (uFz / p.fzmin) : uFz

    V0   = p.longvl
    R0   = p.unloaded_radius
    Fz0  = p.fnomin
    pi0  = p.nompres
    LMY  = p.lmy
    VXLOW = p.vxlow

    QSY1 = p.qsy1; QSY2 = p.qsy2; QSY3 = p.qsy3; QSY4 = p.qsy4
    QSY5 = p.qsy5; QSY6 = p.qsy6; QSY7 = p.qsy7; QSY8 = p.qsy8

    fittyp = p.metadata.fittyp
    is_mf52 = (fittyp == 6 || fittyp == 21)

    My = if is_mf52
        # MF5.2 equation from the MF-Tyre manual
        -R0 * Fz_eff * LMY *
        (QSY1 + QSY2 * (Fx / Fz0) +
         QSY3 * abs(Vcx / V0) + QSY4 * (Vcx / V0)^4)
    else
        # MF6.1 / MF6.2 — Paper Eqn (A48)
        -R0 * Fz0 * LMY *
        (QSY1 + QSY2 * (Fx / Fz0) +
         QSY3 * abs(Vcx / V0) + QSY4 * (Vcx / V0)^4 +
         (QSY5 + QSY6 * (Fz_eff / Fz0)) * gamma^2) *
        ((Fz_eff / Fz0)^QSY7 * (press / pi0)^QSY8)
    end

    # Zero-speed zeroing (must come before low-speed interpolation to avoid ÷0)
    if Vcx == 0.0; return 0.0; end

    # Backward-speed sign flip
    if Vcx < 0.0; My = -My; end

    # Low-speed interpolation (empirically discovered in reference code)
    # The transition region is defined by kappa ∈ [lowLimit, highLimit]
    highLimit = VXLOW / abs(Vcx) - 1.0
    lowLimit  = -1.0 - VXLOW - highLimit

    if kappa >= lowLimit && kappa <= highLimit
        x  = kappa
        x1 = highLimit
        x0 = -1.0
        # Linear interpolation between 0 and π/2, then take sin
        reduction = (0.0 * (x1 - x) + (π / 2.0) * (x - x0)) / (x1 - x0)
        My = My * sin(reduction)
    elseif kappa < lowLimit
        My = -My
    end

    return My
end