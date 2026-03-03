# ==============================================================================
# parse_inputs.jl
#
# Validates, clamps, and pre-processes a single MFInputs point into
# PostProInputs + InternalParams, mirroring Solver.parseInputs in the
# reference MATLAB code.
#
# All output structs are fully concrete (Float64 only) so the hot path
# sees no heap allocations.
# ==============================================================================

# ------------------------------------------------------------------------------
# Output types
# ------------------------------------------------------------------------------

"""
    PostProInputs

Post-processed scalar inputs ready for the force/moment solver.

`u`-prefixed fields are *unlimited* (pre-clamp) copies kept for output
and for quantities (My, Re) that must use raw user values.
"""
struct PostProInputs
    # ── clamped (used in force equations) ──────────────────────────────────
    Fz          ::Float64   # vertical load (≥ 0, saturated to [FZMIN,FZMAX]) [N]
    kappa       ::Float64   # longitudinal slip (saturated)                    [-]
    alpha       ::Float64   # side-slip angle (saturated, low-speed reduced)   [rad]
    gamma       ::Float64   # inclination angle (saturated)                    [rad]
    phit        ::Float64   # turn slip (low-speed reduced)                    [1/m]
    p           ::Float64   # inflation pressure (saturated)                   [Pa]
    # ── unlimited copies ───────────────────────────────────────────────────
    uFz         ::Float64   # raw Fz (output col 3)
    ukappa      ::Float64   # raw kappa, no low-speed reduction (used in My)
    ukappaLow   ::Float64   # kappa with low-speed reduction (output col 7)
    ualpha      ::Float64   # raw alpha (output col 8)
    ugamma      ::Float64   # raw gamma (output col 9)
    uphit       ::Float64   # raw phit  (output col 10)
    uVcx        ::Float64   # raw Vx    (output col 11)
    # ── derived ────────────────────────────────────────────────────────────
    omega       ::Float64   # wheel angular speed                              [rad/s]
    Vsx         ::Float64   # longitudinal slip velocity at S                  [m/s]
    phi         ::Float64   # spin slip (turn-slip mode only)                  [1/m]
    Fz_lowLimit ::Float64   # Fz clamped to [FZMIN, +∞) (used in Mz)         [N]
end

"""
    InternalParams

Solver constants and low-speed reduction scalars.
All epsilon values follow Pacejka Eqn (4.E6a): ε = 1e-6.
"""
struct InternalParams
    epsilonx ::Float64   # protects Bx denominator
    epsilonk ::Float64   # protects Kya denominator
    epsilony ::Float64   # protects By denominator
    epsilonr ::Float64   # (reserved)
    epsilonv ::Float64   # protects Vc denominator
    epsilong ::Float64   # camber-reduction factor (turn-slip only)
    # low-speed reduction scalars (= 1.0 when not in low-speed regime)
    reductionSmooth ::Float64   # smooth cosine ramp  [0,1]
    reductionSharp  ::Float64   # sin(linear*π/2) ramp [0,1]
    reductionLinear ::Float64   # |Vcx|/VXLOW         [0,1]
    reductionLinear_alpha ::Float64   # speedSum/VXLOW for alpha
    # turn-slip zeta2 (scalar; initialised to 1.0, overwritten by Fy0)
    zeta2 ::Float64
end

# ------------------------------------------------------------------------------
# Internal: effective rolling radius used to bootstrap omega
# ------------------------------------------------------------------------------

# Compute Re, Romega, and omega for a single point.
# Called once during parse_inputs when omega is not provided by the user.
# Returns (Re, Romega, omega).
@inline function _compute_re(p::TireParams, uFz::Float64, uVcx::Float64,
                              ukappa::Float64, dpi::Float64,
                              omega_in::Float64) ::Tuple{Float64,Float64,Float64}
    R0   = p.unloaded_radius
    V0   = p.longvl
    Fz0  = p.fnomin
    Cz0  = p.vertical_stiffness
    BREFF = p.breff
    DREFF = p.dreff
    FREFF = p.freff
    Q_RE0 = p.q_re0
    Q_V1  = p.q_v1

    Cz = Cz0 * (1.0 + p.pfz1 * dpi)   # pressure-corrected vertical stiffness

    if omega_in != 0.0
        # omega supplied by user
        omega  = omega_in + 1e-15       # avoid exact zero → Romega = 0
        Romega = R0 * (Q_RE0 + Q_V1 * (omega * R0 / V0)^2)
        Re     = Romega - (Fz0 / Cz) * (DREFF * atan(BREFF * uFz / Fz0) + FREFF * uFz / Fz0)
        return (Re, Romega, omega)
    end

    # Iterative estimation: converges in ≤ 6 iterations for typical inputs
    Re = R0 * 0.965
    Romega = R0 * Q_RE0
    omega  = 0.0
    for _ in 1:20
        Re_prev = Re
        omega   = (ukappa * uVcx + uVcx) / Re     # Eqn (2.5) p.65 Book
        Romega  = R0 * (Q_RE0 + Q_V1 * (omega * R0 / V0)^2)
        Re      = Romega - (Fz0 / Cz) * (DREFF * atan(BREFF * uFz / Fz0) + FREFF * uFz / Fz0)
        abs(Re - Re_prev) < 1e-9 && break
    end
    return (Re, Romega, omega)
end

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

"""
    parse_inputs(p, inp, modes) → (PostProInputs, InternalParams)

Pre-process a single `MFInputs` point: apply physical limits, build the
low-speed reduction factors, compute omega/Re, and pack all derived
quantities needed by the force/moment modules.

Mirrors `Solver.parseInputs` from the reference MATLAB code, steady-state
path only (`userDynamics = 0`).
"""
@inline function parse_inputs(p      ::TireParams,
                               inp    ::MFInputs,
                               modes  ::MFModes) ::Tuple{PostProInputs, InternalParams}

    ε = 1e-6   # universal small factor, Eqn (4.E6a)

    # ── Unpack raw inputs ────────────────────────────────────────────────────
    uFz    = max(inp.Fz, 0.0)   # negative Fz → 0
    kappa  = inp.kappa
    alpha  = inp.alpha
    gamma  = inp.gamma
    phit   = inp.phit
    Vcx    = inp.Vx
    p_in   = inp.pressure
    omega_user = inp.omega

    # Keep unlimited copies for outputs
    ukappaLow = kappa
    ualpha    = alpha
    ugamma    = gamma
    uphit     = phit
    uVcx      = Vcx

    # Zero speed → zero alpha (empirically discovered in reference code)
    if Vcx == 0.0; ualpha = 0.0; end

    # Inflation pressure: use TIR nominal when user supplies 0
    press = (p_in != 0.0) ? p_in : p.inflpres

    # Normalised pressure increment [Eqn (4.E2b)]
    dpi = (press - p.nompres) / p.nompres

    # ── Low-speed reduction factors ──────────────────────────────────────────
    reductionSmooth       = 1.0
    reductionSharp        = 1.0
    reductionLinear       = 1.0
    reductionLinear_alpha = 1.0
    isLowSpeed            = false

    Fz         = uFz
    Fz_lowLimit = uFz

    if modes.use_limits_check
        VXLOW  = p.vxlow
        isLowSpeed = abs(Vcx) <= VXLOW

        if isLowSpeed
            linRed = abs(Vcx / VXLOW)
            reductionLinear = linRed
            Wvlow  = 0.5 * (1.0 + cos(π * (Vcx / VXLOW)))
            reductionSmooth = 1.0 - Wvlow
            reductionSharp  = sin(linRed * (π / 2.0))

            # Low-speed kappa and phit reductions
            kappa       = kappa * linRed
            ukappaLow   = ukappaLow * linRed
            phit        = phit * linRed
        end

        # Negative speed → flip turn slip sign
        if Vcx < 0.0; phit = -phit; end

        # alpha low-speed reduction using combined speed
        Vsy = tan(alpha) * Vcx
        speedSum = abs(Vcx) + abs(Vsy)
        if speedSum < VXLOW
            reductionLinear_alpha = speedSum / VXLOW
            alpha = alpha * reductionLinear_alpha
        end

        # Turn-slip modification (empirically discovered)
        phit = phit * cos(alpha)

        # Saturate to TIR limits
        alpha  = clamp(alpha,  p.alpmin, p.alpmax)
        gamma  = clamp(gamma,  p.cammin, p.cammax)
        Fz     = clamp(uFz,    0.0,      p.fzmax)
        press  = clamp(press,  p.presmin, p.presmax)
        kappa  = clamp(kappa,  p.kpumin, p.kpumax)

        # Fz_lowLimit: Fz clamped to [FZMIN, ∞) — used only in Mz equations
        Fz_lowLimit = (Fz < p.fzmin) ? p.fzmin : Fz
        if uFz <= 0.0; Fz_lowLimit = 0.0; end

    else
        # No limits: set all reductions to 1
        isLowSpeed            = false
        reductionSmooth       = 1.0
        reductionSharp        = 1.0
        reductionLinear       = 1.0
        reductionLinear_alpha = 1.0
    end

    # ── Compute omega and Re ─────────────────────────────────────────────────
    Re, _, omega = _compute_re(p, uFz, Vcx, inp.kappa, dpi, omega_user)

    # Longitudinal slip velocity [Eqn (2.3) p.64 Book]
    Vsx = Vcx - Re * omega

    # ── Turn-slip spin [Eqn (4.76) p.184 Book] ──────────────────────────────
    phi    = 0.0
    epsilong = 0.0

    if modes.use_turn_slip
        Fz0_prime = p.lfzo * p.fnomin
        dfz       = (Fz - Fz0_prime) / Fz0_prime
        epsilong  = p.pecp1 * (1.0 + p.pecp2 * dfz)   # Eqn (4.90)

        Vsy_ts = -tan(alpha) * abs(Vcx)
        V_ts   = sqrt(Vcx^2 + Vsy_ts^2)

        # Speed limit to avoid singularity
        Vc_prime = if abs(Vcx) < p.vxlow
            signVcx = Vcx >= 0.0 ? 1.0 : -1.0
            p.vxlow * signVcx
        else
            V_ts
        end

        psi_dot = -phit * Vc_prime   # Eqn (4.75) rearranged
        phi     = (psi_dot - (1.0 - epsilong) * omega * sin(gamma)) / Vc_prime  # Eqn (4.76)
    end

    pp = PostProInputs(
        Fz, kappa, alpha, gamma, phit, press,
        uFz, inp.kappa, ukappaLow, ualpha, ugamma, uphit, uVcx,
        omega, Vsx, phi, Fz_lowLimit,
    )

    ip = InternalParams(
        ε, ε, ε, ε, ε,     # epsilonx/k/y/r/v
        epsilong,
        reductionSmooth, reductionSharp, reductionLinear, reductionLinear_alpha,
        1.0,               # zeta2 — overwritten by Fy0
    )

    return (pp, ip)
end