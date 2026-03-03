# ==============================================================================
# geometry.jl
#
# Tyre geometry calculations:
#   • calc_Re    — effective rolling radius, Romega, omega
#   • calc_rho_Rl — tyre deflection, loaded radius, vertical stiffness
#   • calc_contact_patch — contact half-lengths a, b
#
# MF5.2/6.1 uses the quadratic Eqn (A3.3); MF6.2 solves iteratively via
# the secant method since Rl appears implicitly in Fz(Rl).
#
# References:
#   Pacejka (2012) Book App. 3, Eqns (A3.2–A3.3)
#   Besselink et al. paper, Eqns (1,4–10)
#   MF-Swift 6.2 equation manual (2013-07-06)
# ==============================================================================

# ------------------------------------------------------------------------------
# Effective rolling radius
# (Re calculation is duplicated from parse_inputs._compute_re to expose
#  all three outputs; the per-evaluation call inside parse_inputs avoids
#  paying this cost twice for the common steady-state path.)
# ------------------------------------------------------------------------------

"""
    calc_Re(p, pp, dpi) → (Re, Romega, omega)

Effective rolling radius (m), free-spinning radius (m), and angular speed
(rad/s).  When `pp.omega ≠ 0`, omega is taken directly; otherwise it is
estimated iteratively from Vcx and kappa.
"""
@inline @fastmath function calc_Re(p   ::TireParams,
                                    pp  ::PostProInputs,
                                    dpi ::Float64) ::Tuple{Float64,Float64,Float64}

    Vcx    = pp.uVcx
    uFz    = pp.uFz
    ukappa = pp.ukappa   # truly unlimited kappa (used in omega estimation, matches reference)

    R0   = p.unloaded_radius
    V0   = p.longvl
    Fz0  = p.fnomin
    Cz0  = p.vertical_stiffness
    BREFF = p.breff
    DREFF = p.dreff
    FREFF = p.freff
    Q_RE0 = p.q_re0
    Q_V1  = p.q_v1

    Cz = Cz0 * (1.0 + p.pfz1 * dpi)   # Eqn (5) paper — pressure-corrected Cz

    if pp.omega != 0.0
        omega  = pp.omega + 1e-15
        Romega = R0 * (Q_RE0 + Q_V1 * (omega * R0 / V0)^2)            # Eqn (1)
        Re     = Romega - (Fz0 / Cz) *
                 (DREFF * atan(BREFF * uFz / Fz0) + FREFF * uFz / Fz0)  # Eqn (7)
        return (Re, Romega, omega)
    end

    Re     = R0 * 0.965   # initial guess
    Romega = R0 * Q_RE0
    omega  = 0.0
    for _ in 1:20
        Re_prev = Re
        omega   = (ukappa * Vcx + Vcx) / Re                            # Eqn (2.5)
        Romega  = R0 * (Q_RE0 + Q_V1 * (omega * R0 / V0)^2)
        Re      = Romega - (Fz0 / Cz) *
                  (DREFF * atan(BREFF * uFz / Fz0) + FREFF * uFz / Fz0)
        abs(Re - Re_prev) < 1e-9 && break
    end
    return (Re, Romega, omega)
end

# ------------------------------------------------------------------------------
# MF6.2 — vertical force from loaded radius (needed for secant solver)
# ------------------------------------------------------------------------------

@inline @fastmath function _fz62(p      ::TireParams,
                                  gamma  ::Float64,
                                  omega  ::Float64,
                                  Romega ::Float64,
                                  dpi    ::Float64,
                                  Rl     ::Float64,
                                  Fx     ::Float64,
                                  Fy     ::Float64) ::Tuple{Float64,Float64}

    V0   = p.longvl
    R0   = p.unloaded_radius
    AR   = p.aspect_ratio
    W    = p.width
    Fz0  = p.fnomin
    Cz0  = p.vertical_stiffness
    Q_V2   = p.q_v2;   Q_FZ2  = p.q_fz2
    Q_FCX  = p.q_fcx;  Q_FCY  = p.q_fcy
    PFZ1   = p.pfz1
    Q_FCY2 = p.q_fcy2
    Q_CAM1 = p.q_cam1; Q_CAM2 = p.q_cam2; Q_CAM3 = p.q_cam3
    Q_FYS1 = p.q_fys1; Q_FYS2 = p.q_fys2; Q_FYS3 = p.q_fys3

    Q_FZ1 = sqrt(max((Cz0 * R0 / Fz0)^2 - 4.0 * Q_FZ2, 0.0))  # rearranged Eqn (4)

    # Asymmetric camber–lateral-force effect
    Sfyg = (Q_FYS1 + Q_FYS2 * (Rl / Romega) + Q_FYS3 * (Rl / Romega)^2) * gamma

    # Free-rolling deflection
    rho_zfr = max(Romega - Rl, 0.0)

    # Reference tread width
    rtw = (1.075 - 0.5 * AR) * W

    # Camber deflection — guard against γ=0 → 0/0
    cam_num = (Q_CAM1 * Rl    + Q_CAM2 * Rl^2)    * gamma
    cam_den = (Q_CAM1 * Romega + Q_CAM2 * Romega^2) * gamma
    rho_zg  = if abs(cam_den) > 1e-15
        (cam_num / cam_den)^2 * (rtw / 8.0) * abs(tan(gamma)) -
        Q_CAM3 * rho_zfr * abs(gamma)
    else
        0.0
    end

    rho_z = max(rho_zfr + rho_zg, 1e-12)

    # Correction factor
    fcorr = (1.0 + Q_V2 * (R0 / V0) * abs(omega) -
             (Q_FCX * Fx / Fz0)^2 -
             (rho_z / R0)^Q_FCY2 * (Q_FCY * (Fy - Sfyg) / Fz0)^2) *
            (1.0 + PFZ1 * dpi)

    Fz_calc = fcorr * (Q_FZ1 * (rho_z / R0) + Q_FZ2 * (rho_z / R0)^2) * Fz0

    return (Fz_calc, rho_z)
end

# ------------------------------------------------------------------------------
# Secant solver — scalar, in-place, for MF6.2 Rl
# ------------------------------------------------------------------------------

@inline function _secant_Rl(p      ::TireParams,
                              gamma  ::Float64,
                              omega  ::Float64,
                              Romega ::Float64,
                              dpi    ::Float64,
                              Fx     ::Float64,
                              Fy     ::Float64,
                              Fz_tgt ::Float64) ::Float64

    R0  = p.unloaded_radius
    tol = p.fnomin * 1e-6

    x0 = R0 * 0.95
    x1 = R0

    fz0, _ = _fz62(p, gamma, omega, Romega, dpi, x0, Fx, Fy)
    fz1, _ = _fz62(p, gamma, omega, Romega, dpi, x1, Fx, Fy)
    y0 = fz0 - Fz_tgt
    y1 = fz1 - Fz_tgt

    x = x1
    for _ in 3:20
        dy = y1 - y0
        abs(dy) < 1e-15 && break
        try_x = x1 - (x1 - x0) / dy * y1
        x  = try_x
        fz, _ = _fz62(p, gamma, omega, Romega, dpi, x, Fx, Fy)
        y  = fz - Fz_tgt
        abs(y) <= tol && break
        if abs(y1) < abs(y0)
            x0 = x1; y0 = y1
        end
        x1 = x; y1 = y
    end
    return x
end

# ------------------------------------------------------------------------------
# rho, Rl, Cz — MF5.2 / MF6.1
# ------------------------------------------------------------------------------

@inline @fastmath function _rho_Rl_61(p      ::TireParams,
                                       pp     ::PostProInputs,
                                       dpi    ::Float64,
                                       omega  ::Float64,
                                       Romega ::Float64,
                                       Fx     ::Float64,
                                       Fy     ::Float64) ::Tuple{Float64,Float64,Float64}

    Fz  = pp.Fz_lowLimit   # uses low-limit Fz as in reference code
    R0  = p.unloaded_radius
    V0  = p.longvl
    Fz0 = p.fnomin
    Cz0 = p.vertical_stiffness
    Q_V2  = p.q_v2;  Q_FZ2 = p.q_fz2
    Q_FCX = p.q_fcx; Q_FCY = p.q_fcy
    PFZ1  = p.pfz1

    Q_FZ1 = sqrt(max((Cz0 * R0 / Fz0)^2 - 4.0 * Q_FZ2, 0.0))

    ext = (1.0 + Q_V2 * (R0 / V0) * abs(omega) -
           (Q_FCX * Fx / Fz0)^2 -
           (Q_FCY * Fy / Fz0)^2) * (1.0 + PFZ1 * dpi) * Fz0

    # Quadratic: A·x² + B·x + C = 0  where x = rho/R0
    A =  Q_FZ2
    B =  Q_FZ1
    C = -(Fz / max(ext, 1e-15))

    disc = B^2 - 4.0 * A * C
    x    = (-B + sqrt(max(disc, 0.0))) / (2.0 * A)

    rho = max(x * R0, 0.0)
    Rl  = max(Romega - rho, 0.0)
    rho = (Fz == 0.0) ? 1e-6 : rho

    Cz  = pp.Fz_lowLimit / rho

    return (rho, Rl, Cz)
end

# ------------------------------------------------------------------------------
# rho, Rl, Cz — MF6.2
# ------------------------------------------------------------------------------

@inline function _rho_Rl_62(p      ::TireParams,
                              pp     ::PostProInputs,
                              dpi    ::Float64,
                              omega  ::Float64,
                              Romega ::Float64,
                              Fx     ::Float64,
                              Fy     ::Float64) ::Tuple{Float64,Float64,Float64}

    uFz   = pp.uFz
    gamma = pp.gamma
    Fz_tgt = uFz   # MF6.2 uses uFz as target

    Rl = _secant_Rl(p, gamma, omega, Romega, dpi, Fx, Fy, Fz_tgt)

    _, rho = _fz62(p, gamma, omega, Romega, dpi, Rl, Fx, Fy)

    # Guard zero-force case
    rho = (pp.uFz == 0.0) ? 1e-6 : rho
    Cz  = pp.uFz / rho

    return (rho, Rl, Cz)
end

# ------------------------------------------------------------------------------
# Public dispatcher
# ------------------------------------------------------------------------------

"""
    calc_rho_Rl(p, pp, dpi, omega, Romega, Fx, Fy) → (rho, Rl, Cz)

Tyre deflection (m), loaded radius (m), and vertical stiffness (N/m).
Dispatches to MF6.2 secant solver or the MF5.2/6.1 analytic solution.
"""
@inline function calc_rho_Rl(p      ::TireParams,
                               pp     ::PostProInputs,
                               dpi    ::Float64,
                               omega  ::Float64,
                               Romega ::Float64,
                               Fx     ::Float64,
                               Fy     ::Float64) ::Tuple{Float64,Float64,Float64}
    if p.metadata.fittyp == 62
        return _rho_Rl_62(p, pp, dpi, omega, Romega, Fx, Fy)
    else
        return _rho_Rl_61(p, pp, dpi, omega, Romega, Fx, Fy)
    end
end

# ------------------------------------------------------------------------------
# Contact patch half-lengths a and b
# ------------------------------------------------------------------------------

"""
    calc_contact_patch(p, pp, dpi) → (a, b, NCz)

Contact patch semi-length (m), semi-width (m), and pressure-corrected
vertical stiffness NCz (N/m).

MF5.2 uses Q_A1/Q_A2 (with default derivation when both are zero).
MF6.1/6.2 use Q_RA1/Q_RA2 and Q_RB1/Q_RB2.
"""
@inline @fastmath function calc_contact_patch(p   ::TireParams,
                                               pp  ::PostProInputs,
                                               dpi ::Float64) ::Tuple{Float64,Float64,Float64}

    uFz = pp.uFz
    R0  = p.unloaded_radius
    w   = p.width
    Cz0 = p.vertical_stiffness
    Fz0 = p.fnomin
    PFZ1 = p.pfz1

    NCz = Cz0 * (1.0 + PFZ1 * dpi)   # Eqn (5) paper

    # Bottoming check
    Rrim  = p.rim_radius
    Dbtm  = p.bottom_offst
    Rl_approx = R0 - uFz / NCz
    Fz_eff = if Rl_approx - (Rrim + Dbtm) < 0.0
        (R0 - Rrim - Dbtm) * NCz
    else
        uFz
    end

    fittyp = p.metadata.fittyp
    if fittyp == 6 || fittyp == 21
        # MF5.2
        Q_A1 = p.q_a1; Q_A2 = p.q_a2
        if Q_A1 == 0.0 && Q_A2 == 0.0
            y    = log10(R0 * Cz0 / Fz0)
            Q_A1 = -0.0388 * y^3 + 0.2509 * y^2 - 0.6283 * y + 0.6279
            Q_A2 = 1.693 * Q_A1^2
        end
        a = R0 * (Q_A2 * (Fz_eff / Fz0) + Q_A1 * sqrt(Fz_eff / Fz0))
        b = w / 2.0
    else
        # MF6.1 / MF6.2 — Eqns (9) and (10) of paper
        Q_RA1 = p.q_ra1; Q_RA2 = p.q_ra2
        Q_RB1 = p.q_rb1; Q_RB2 = p.q_rb2
        nratio = Fz_eff / (NCz * R0)
        a = R0 * (Q_RA2 * nratio + Q_RA1 * sqrt(nratio))
        b = w  * (Q_RB2 * nratio + Q_RB1 * nratio^(1.0/3.0))
    end

    return (a, b, NCz)
end