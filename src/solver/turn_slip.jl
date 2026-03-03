# ==============================================================================
# turn_slip.jl
#
# Turn-slip (spin-slip) ζ factors.
# These are called from Fx0/Fy0/Mz0 when modes.use_turn_slip == true.
#
# All equations reference Pacejka 2012 book, Chapter 4.3.3 (pp. 183-190).
# ==============================================================================

"""
    TurnSlipFactors

The ζ₀–ζ₈ weighting factors for turn-slip influence on forces and moments.
When `use_turn_slip = false`, every ζ is exactly 1 (except ζ₀ = 0 per book).
"""
struct TurnSlipFactors
    zeta0 ::Float64   # horizontal-shift factor for Fy: 0 (turn slip) / 1 (no)
    zeta1 ::Float64   # peak-Fx reduction due to spin        Eqn (4.105)
    zeta2 ::Float64   # peak-Fy reduction due to spin        Eqn (4.77)
    zeta3 ::Float64   # cornering-stiffness reduction        Eqn (4.79)
    zeta4 ::Float64   # Fy horizontal-shift factor           Eqn (4.84)
    zeta5 ::Float64   # pneumatic-trail reduction            Eqn (4.91)
    zeta6 ::Float64   # residual-torque reduction            Eqn (4.102)
    zeta7 ::Float64   # Mzr reduction (acos based)           manual
    zeta8 ::Float64   # Mzr additive offset                  manual
end

# Convenience: unit turn-slip factors (no turn slip)
@inline _no_turn_slip() = TurnSlipFactors(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

# ------------------------------------------------------------------------------
# zeta1 — called from Fx0
# ------------------------------------------------------------------------------

"""
    calc_zeta1(p, pp, dfz) → Float64

ζ₁ = cos(atan(Bxp · R₀ · φ))   [Eqn (4.105)]

Only used when `use_turn_slip = true`.
"""
@inline @fastmath function calc_zeta1(p   ::TireParams,
                                       pp  ::PostProInputs,
                                       dfz ::Float64) ::Float64
    R0    = p.unloaded_radius
    PDXP1 = p.pdxp1
    PDXP2 = p.pdxp2
    PDXP3 = p.pdxp3

    Bxp = PDXP1 * (1.0 + PDXP2 * dfz) * cos(atan(PDXP3 * pp.kappa))  # Eqn (4.106)
    return cos(atan(Bxp * R0 * pp.phi))                                 # Eqn (4.105)
end

# ------------------------------------------------------------------------------
# zeta2 & zeta3 — called from Fy0
# ------------------------------------------------------------------------------

"""
    calc_zeta2_zeta3(p, pp, dfz) → (Float64, Float64)

ζ₂ = cos(atan(Byp·(R₀|φ| + PDYP4·√(R₀|φ|))))   [Eqn (4.77)]
ζ₃ = cos(atan(PKYP1·R₀²·φ²))                     [Eqn (4.79)]
"""
@inline @fastmath function calc_zeta2_zeta3(p   ::TireParams,
                                              pp  ::PostProInputs,
                                              dfz ::Float64) ::Tuple{Float64,Float64}
    R0    = p.unloaded_radius
    PKYP1 = p.pkyp1
    PDYP1 = p.pdyp1
    PDYP2 = p.pdyp2
    PDYP3 = p.pdyp3
    PDYP4 = p.pdyp4

    phi = pp.phi

    zeta3 = cos(atan(PKYP1 * R0^2 * phi^2))                                   # Eqn (4.79)

    Byp   = PDYP1 * (1.0 + PDYP2 * dfz) * cos(atan(PDYP3 * tan(pp.alpha)))   # Eqn (4.78)
    R0phi = R0 * abs(phi)
    zeta2 = cos(atan(Byp * (R0phi + PDYP4 * sqrt(R0phi))))                    # Eqn (4.77)

    return (zeta2, zeta3)
end

# ------------------------------------------------------------------------------
# Mz0 turn-slip terms — called from Mz0
# ------------------------------------------------------------------------------

"""
    calc_mz0_turnslip(p, pp, ip, modes, iv, Kya, Kyg0, muy, Fy0, Gyk, Fz0_prime)
        → (zeta5, zeta6, zeta7, zeta8)

Computes the four ζ factors needed for the Mz0 aligning torque under turn slip.

Arguments `Fy0`, `Gyk` are already computed by Fy0/combined-slip at this call site.
"""
@inline @fastmath function calc_mz0_turnslip(p        ::TireParams,
                                               pp       ::PostProInputs,
                                               ip       ::InternalParams,
                                               iv       ::IncrVars,
                                               Kya      ::Float64,
                                               Kyg0     ::Float64,
                                               muy      ::Float64,
                                               Fy0      ::Float64,
                                               Gyk      ::Float64,
                                               Fz0_prime::Float64) ::Tuple{Float64,Float64,Float64,Float64}

    R0    = p.unloaded_radius
    QDTP1 = p.qdtp1
    QBRP1 = p.qbrp1
    QCRP1 = p.qcrp1
    QDRP1 = p.qdrp1
    QCRP2 = p.qcrp2
    QDZ8  = p.qdz8
    QDZ9  = p.qdz9
    QDZ10 = p.qdz10
    QDZ11 = p.qdz11
    LMP   = p.lmp
    LKZC  = p.lkzc

    Fz        = pp.Fz
    phi       = pp.phi
    phit      = pp.phit
    gamma     = pp.gamma
    epsilong  = ip.epsilong
    dfz       = iv.dfz
    dpi       = iv.dpi

    zeta5 = cos(atan(QDTP1 * R0 * phi))   # Eqn (4.91)
    zeta6 = cos(atan(QBRP1 * R0 * phi))   # Eqn (4.102)

    # Mzp_inf must be > 0  [Eqn (4.95)]
    Mzp_inf = max(QCRP1 * abs(muy) * R0 * Fz * sqrt(Fz / Fz0_prime) * LMP, 1e-6)

    CDrp = QDRP1                                    # Eqn (4.96) — (>0)
    DDrp = Mzp_inf / sin(0.5 * π * CDrp)           # Eqn (4.94)

    Kzgr0 = Fz * R0 * (QDZ8 * QDZ9 * dfz +
             (QDZ10 + QDZ11 * dfz * abs(gamma))) * LKZC   # Eqn (4.99)

    BDrp  = Kzgr0 / (CDrp * DDrp * (1.0 - epsilong))     # manual
    Drp   = DDrp * sin(CDrp * atan(BDrp * R0 * phit))     # manual

    Mzp90 = Mzp_inf * (2.0 / π) * atan(QCRP2 * R0 * abs(phit)) * Gyk   # Eqn (4.103)

    # Guard: acos argument must be in [-1, 1]
    acos_arg = clamp(Mzp90 / max(abs(Drp), 1e-15), -1.0, 1.0)
    zeta7 = (2.0 / π) * acos(acos_arg)   # manual
    zeta8 = 1.0 + Drp                    # manual

    return (zeta5, zeta6, zeta7, zeta8)
end