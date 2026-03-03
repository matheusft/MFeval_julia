# ==============================================================================
# basic_vars.jl
#
# Shared pre-computations performed once per evaluation before any
# force/moment module is called.  Mirrors Solver.calculateBasic.
#
# All inputs and outputs are plain Float64 scalars — no allocations.
# ==============================================================================

"""
    StarVars

Starred quantities (alpha_star, gamma_star) and friction scaling factors
used throughout the force/moment equations.
"""
struct StarVars
    alpha_star  ::Float64   # side-slip measure: tan(α)·sign(Vcx) or α   [rad or -]
    gamma_star  ::Float64   # camber measure: sin(γ) or γ                 [rad or -]
    LMUX_star   ::Float64   # friction scale with speed decay (Eqn 4.E7)  [-]
    LMUY_star   ::Float64   # friction scale with speed decay (Eqn 4.E7)  [-]
end

"""
    PrimeVars

Primed quantities used by the force/moment modules.
"""
struct PrimeVars
    Fz0_prime   ::Float64   # LFZO·Fz0 (Eqn 4.E1)                       [N]
    alpha_prime ::Float64   # acos(Vcx/Vc') (Eqn 4.E6)                   [rad]
    LMUX_prime  ::Float64   # digressive friction factor X (Eqn 4.E8)     [-]
    LMUY_prime  ::Float64   # digressive friction factor Y (Eqn 4.E8)     [-]
end

"""
    IncrVars

Normalised increments reused everywhere.
"""
struct IncrVars
    dfz ::Float64   # (Fz − Fz0')/Fz0'  (Eqn 4.E2a)  [-]
    dpi ::Float64   # (p  − π0)/π0       (Eqn 4.E2b)  [-]
end

"""
    SlipVelocities

Slip-point velocities (used only in turn-slip and relaxation modules).
"""
struct SlipVelocities
    Vsx ::Float64   # longitudinal slip velocity at S (Eqn 4.E5)  [m/s]
    Vsy ::Float64   # lateral    slip velocity at S (Eqn 2.12)    [m/s]
end

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

"""
    calc_basic_vars(p, pp, modes) → (StarVars, PrimeVars, IncrVars, SlipVelocities)

Compute all shared pre-variables needed before any force/moment calculation.
Mirrors `Solver.calculateBasic` from the reference MATLAB code.
"""
@inline @fastmath function calc_basic_vars(p     ::TireParams,
                                            pp    ::PostProInputs,
                                            ip    ::InternalParams,
                                            modes ::MFModes) ::Tuple{StarVars, PrimeVars, IncrVars, SlipVelocities}

    # ── Unpack frequently used quantities ────────────────────────────────────
    V0   = p.longvl
    Fz0  = p.fnomin
    LFZO = p.lfzo
    LMUX = p.lmux
    LMUY = p.lmuy

    Fz  = pp.Fz
    Vcx = pp.uVcx

    # Constant from Pacejka 2012 — Amu = 1 gives perfect match with TNO
    Amu = 1.0

    # ── Slip velocities ──────────────────────────────────────────────────────
    # Vsx = −κ·|Vcx|                      Eqn (4.E5)
    # Vsy =  tan(α)·|Vcx|                 Eqn (2.12) / (4.E3)
    Vsx = -pp.kappa * abs(Vcx)
    Vsy =  tan(pp.alpha) * abs(Vcx)

    # Total slip speed at contact point S  Eqn (3.39)
    Vs = sqrt(Vsx^2 + Vsy^2)

    # Contact-centre velocity
    Vcy = Vsy    # assumption: p.67 Book, above Eqn (2.11)
    Vc  = sqrt(Vcx^2 + Vcy^2)

    # ── Fz0_prime and increments ─────────────────────────────────────────────
    Fz0_prime = LFZO * Fz0                           # Eqn (4.E1)
    dfz       = (Fz - Fz0_prime) / Fz0_prime         # Eqn (4.E2a)
    dpi       = (pp.p - p.nompres) / p.nompres        # Eqn (4.E2b)

    # ── alpha_star / gamma_star ───────────────────────────────────────────────
    if modes.use_alpha_star
        alpha_star = tan(pp.alpha) * (Vcx >= 0.0 ? 1.0 : -1.0)   # Eqn (4.E3)
        gamma_star = sin(pp.gamma)                                  # Eqn (4.E4)
    else
        alpha_star = pp.alpha
        gamma_star = pp.gamma
    end

    # ── alpha_prime (for aligning torque at high slip) ───────────────────────
    # Vc_prime = Vc + ε·sign(Vc)  Eqn (4.E6a)
    signVc   = Vc >= 0.0 ? 1.0 : -1.0
    Vc_prime = Vc + ip.epsilonv * signVc
    alpha_prime = acos(clamp(Vcx / Vc_prime, -1.0, 1.0))   # Eqn (4.E6)

    # ── Friction scaling ─────────────────────────────────────────────────────
    # Speed-decaying friction (LMUV = 0 by default — not in standard TIR files)
    LMUV = 0.0
    denom = 1.0 + LMUV * Vs / V0
    LMUX_star = LMUX / denom    # Eqn (4.E7)
    LMUY_star = LMUY / denom    # Eqn (4.E7)

    # Digressive friction factor  Eqn (4.E8)
    # Amu = 1 → LMUX_prime = LMUX_star (no digressivity) — matches TNO
    LMUX_prime = Amu * LMUX_star / (1.0 + (Amu - 1.0) * LMUX_star)
    LMUY_prime = Amu * LMUY_star / (1.0 + (Amu - 1.0) * LMUY_star)

    sv = StarVars(alpha_star, gamma_star, LMUX_star, LMUY_star)
    pv = PrimeVars(Fz0_prime, alpha_prime, LMUX_prime, LMUY_prime)
    iv = IncrVars(dfz, dpi)
    sl = SlipVelocities(Vsx, Vsy)

    return (sv, pv, iv, sl)
end