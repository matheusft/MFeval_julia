# ==============================================================================
# Fy0.jl
#
# Pure lateral force Fy0, lateral friction coefficient muy, cornering
# stiffness Kya, and intermediate quantities needed by Mz0.
#
# Mirrors Solver.calculateFy0 from the reference MATLAB code.
#
# Reference: Pacejka (2012) Book, pp. 177-183 (4.E19–4.E30)
#            MF5.2 equations from the MF-Tyre equation manual.
# ==============================================================================

"""
    Fy0Intermediates

Quantities produced by `calc_Fy0` and consumed by `calc_Mz0` / `calc_Mz`.
Bundled to avoid recomputing them.
"""
struct Fy0Intermediates
    SHy  ::Float64   # horizontal shift of Fy0 curve     (4.E27)
    SVy  ::Float64   # vertical shift of Fy0 curve       (4.E29)
    By   ::Float64   # stiffness factor                  (4.E26)
    Cy   ::Float64   # shape factor                      (4.E21)
    Kya  ::Float64   # cornering stiffness               (4.E25)
    Kyg0 ::Float64   # camber stiffness at γ=0           (4.E30 or MF5.2)
    SVyg ::Float64   # camber-induced vertical shift     (4.E28)
end

# ------------------------------------------------------------------------------
# Helper: evaluate Fy0 at gamma=0, phit=0 (needed by Mz0 for Mzo')
# The caller is responsible for passing zeroed gamma_star.
# ------------------------------------------------------------------------------

"""
    _fy0_gamma0(p, pp_g0, ip, modes, sv_g0, pv, iv, dfz, dpi, Fz0_prime, LMUY_star, LMUY_prime)
        → (Fy0_g0, muy_g0)

Evaluate pure lateral force with γ = 0 and no turn-slip.
Used internally by `calc_Fy0` and `calc_Mz0`.
"""
@inline @fastmath function _fy0_gamma0(p         ::TireParams,
                                        alpha_star ::Float64,
                                        Fz         ::Float64,
                                        Vcx        ::Float64,
                                        ip         ::InternalParams,
                                        pv         ::PrimeVars,
                                        iv         ::IncrVars,
                                        modes      ::MFModes) ::Tuple{Float64,Float64}

    dfz        = iv.dfz
    dpi        = iv.dpi
    Fz0_prime  = pv.Fz0_prime
    LMUY_star  = p.lmuy   # LMUV=0 → LMUY_star = LMUY
    LMUY_prime = pv.LMUY_prime
    epsilony   = ip.epsilony
    epsilonk   = ip.epsilonk

    LCY  = p.lcy
    LEY  = p.ley
    LKY  = p.lky
    LHY  = p.lhy
    LVY  = p.lvy
    LKYC = p.lkyc

    PCY1 = p.pcy1
    PDY1 = p.pdy1; PDY2 = p.pdy2; PDY3 = p.pdy3
    PEY1 = p.pey1; PEY2 = p.pey2; PEY3 = p.pey3; PEY4 = p.pey4; PEY5 = p.pey5
    PKY1 = p.pky1; PKY2 = p.pky2; PKY3 = p.pky3; PKY4 = p.pky4; PKY5 = p.pky5
    PHY1 = p.phy1; PHY2 = p.phy2
    PVY1 = p.pvy1; PVY2 = p.pvy2; PVY3 = p.pvy3; PVY4 = p.pvy4
    PPY1 = p.ppy1; PPY2 = p.ppy2; PPY3 = p.ppy3; PPY4 = p.ppy4; PPY5 = p.ppy5

    # γ=0, so gamma_star = 0, SVyg = 0, zeta2 = 1, zeta3 = 1
    Kya = PKY1 * Fz0_prime * (1.0 + PPY1 * dpi) *
          sin(PKY4 * atan((Fz / Fz0_prime) /
          ((PKY2 + PKY5 * 0.0) * (1.0 + PPY2 * dpi)))) * LKY   # (4.E25) at γ=0

    SHy = (PHY1 + PHY2 * dfz) * LHY                             # (4.E27) at γ=0
    SVy = Fz * (PVY1 + PVY2 * dfz) * LVY * LMUY_prime           # (4.E29) at γ=0

    alphay = alpha_star + SHy                                    # (4.E20)
    Cy     = PCY1 * LCY                                          # (4.E21)
    muy    = (PDY1 + PDY2 * dfz) *
             (1.0 + PPY3 * dpi + PPY4 * dpi^2) * LMUY_star      # (4.E23) at γ=0
    muy    = (Fz == 0.0) ? 0.0 : muy
    Dy     = muy * Fz                                            # (4.E22) at γ=0

    Ey     = (PEY1 + PEY2 * dfz) *
             (1.0 - PEY3 * (alphay >= 0.0 ? 1.0 : -1.0)) * p.ley  # (4.E24) at γ=0
    Ey     = min(Ey, 1.0)

    signDy = (Dy >= 0.0) ? 1.0 : -1.0
    By     = Kya / (Cy * Dy + epsilony * signDy)                 # (4.E26)

    Fy0 = Dy * sin(Cy * atan(By * alphay - Ey * (By * alphay - atan(By * alphay)))) + SVy  # (4.E19)

    if modes.use_alpha_star && Vcx < 0.0; Fy0 = -Fy0; end

    return (Fy0, muy)
end

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

"""
    calc_Fy0(p, pp, ip, modes, sv, pv, iv) → (Fy0, muy, Fy0Intermediates)

Pure lateral force (N), lateral friction coefficient (-), and a bundle of
intermediate quantities needed by the aligning moment modules.

Handles MF5.2 vs MF6.1/6.2 branching, turn-slip ζ₂/ζ₃, and low-speed
smoothing.  `ip.zeta2` is **updated in-place** (returned via new struct).
"""
@inline @fastmath function calc_Fy0(p     ::TireParams,
                                     pp    ::PostProInputs,
                                     ip    ::InternalParams,
                                     modes ::MFModes,
                                     sv    ::StarVars,
                                     pv    ::PrimeVars,
                                     iv    ::IncrVars) ::Tuple{Float64,Float64,Fy0Intermediates,InternalParams}

    Fz         = pp.Fz
    Vcx        = pp.uVcx
    dfz        = iv.dfz
    dpi        = iv.dpi
    Fz0_prime  = pv.Fz0_prime
    LMUY_star  = sv.LMUY_star
    alpha_star = sv.alpha_star
    gamma_star = sv.gamma_star
    LMUY_prime = pv.LMUY_prime
    epsilony   = ip.epsilony
    epsilonk   = ip.epsilonk

    # ── Scaling coefficients ─────────────────────────────────────────────────
    LCY  = p.lcy
    LKY  = p.lky
    LHY  = p.lhy
    LVY  = p.lvy
    LKYC = p.lkyc

    # ── Lateral force coefficients ───────────────────────────────────────────
    PCY1 = p.pcy1
    PDY1 = p.pdy1; PDY2 = p.pdy2; PDY3 = p.pdy3
    PEY1 = p.pey1; PEY2 = p.pey2; PEY3 = p.pey3; PEY4 = p.pey4; PEY5 = p.pey5
    PKY1 = p.pky1; PKY2 = p.pky2; PKY3 = p.pky3; PKY4 = p.pky4; PKY5 = p.pky5
    PKY6 = p.pky6; PKY7 = p.pky7
    PHY1 = p.phy1; PHY2 = p.phy2
    PVY1 = p.pvy1; PVY2 = p.pvy2; PVY3 = p.pvy3; PVY4 = p.pvy4
    PPY1 = p.ppy1; PPY2 = p.ppy2; PPY3 = p.ppy3; PPY4 = p.ppy4; PPY5 = p.ppy5

    # ── Turn-slip ζ₂ / ζ₃ ────────────────────────────────────────────────────
    zeta2, zeta3 = modes.use_turn_slip ?
                   calc_zeta2_zeta3(p, pp, dfz) : (1.0, 1.0)
    zeta0 = modes.use_turn_slip ? 0.0 : 1.0

    # ── Cornering stiffness ───────────────────────────────────────────────────
    Kya = PKY1 * Fz0_prime * (1.0 + PPY1 * dpi) *
          (1.0 - PKY3 * abs(gamma_star)) *
          sin(PKY4 * atan((Fz / Fz0_prime) /
          ((PKY2 + PKY5 * gamma_star^2) * (1.0 + PPY2 * dpi)))) *
          zeta3 * LKY                                                    # (4.E25)

    SVyg = Fz * (PVY3 + PVY4 * dfz) * gamma_star * LKYC * LMUY_prime * zeta2  # (4.E28)

    # ── MF5.2 vs MF6.1/6.2 branching ────────────────────────────────────────
    fittyp = p.metadata.fittyp
    is_mf52 = (fittyp == 6 || fittyp == 21)

    Kyg0 = if is_mf52
        PHY3  = p.phy3
        Kya0  = PKY1 * Fz0_prime *
                sin(PKY4 * atan(Fz / (PKY2 * Fz0_prime))) * LKYC   # MF5.2 simplified
        PHY3 * Kya0 + Fz * (PVY3 + PVY4 * dfz) * LKYC             # MF5.2 Kyg0
    else
        Fz * (PKY6 + PKY7 * dfz) * (1.0 + PPY5 * dpi) * LKYC      # (4.E30)
    end

    # ── Turn-slip SHyp (horizontal shift due to spin) ────────────────────────
    SHyp = 0.0
    zeta4 = 1.0

    if modes.use_turn_slip
        PHYP1 = p.phyp1; PHYP2 = p.phyp2; PHYP3 = p.phyp3; PHYP4 = p.phyp4
        R0    = p.unloaded_radius
        epsilong = ip.epsilong
        phi   = pp.phi

        # Kya0: cornering stiffness at γ=0
        Kya0  = PKY1 * Fz0_prime * (1.0 + PPY1 * dpi) *
                sin(PKY4 * atan((Fz / Fz0_prime) /
                ((PKY2) * (1.0 + PPY2 * dpi)))) * zeta3 * LKY

        signKya  = (Kya  >= 0.0) ? 1.0 : -1.0
        signKya0 = (Kya0 >= 0.0) ? 1.0 : -1.0
        Kya_prime  = Kya  + epsilonk * signKya      # (4.E39)
        Kyao_prime = Kya0 + epsilonk * signKya0

        KyRp0  = Kyg0 / (1.0 - epsilong)            # Eqn (4.89)

        CHyp   = PHYP1                               # (4.E85) (>0)
        DHyp   = (PHYP2 + PHYP3 * dfz) * (Vcx >= 0.0 ? 1.0 : -1.0)   # (4.E86)
        EHyp   = min(PHYP4, 1.0)                    # (4.E87) ≤1

        BHyp   = KyRp0 / (CHyp * DHyp * Kyao_prime)  # (4.E88)
        SHyp   = DHyp * sin(CHyp * atan(BHyp * R0 * phi -
                 EHyp * (BHyp * R0 * phi - atan(BHyp * R0 * phi)))) *
                 (Vcx >= 0.0 ? 1.0 : -1.0)            # (4.E80)

        zeta4  = 1.0 + SHyp - SVyg / Kya_prime        # (4.E84)
    end

    # ── Horizontal shift SHy ─────────────────────────────────────────────────
    SHy = if is_mf52
        PHY3 = p.phy3
        (PHY1 + PHY2 * dfz) * LHY + PHY3 * gamma_star * LKYC   # MF5.2 manual
    else
        signKya_s = (Kya >= 0.0) ? 1.0 : -1.0
        (PHY1 + PHY2 * dfz) * LHY +
        ((Kyg0 * gamma_star - SVyg) / (Kya + epsilonk * signKya_s)) *
        zeta0 + zeta4 - 1.0                                     # (4.E27)
    end

    # ── Vertical shift SVy ────────────────────────────────────────────────────
    SVy = Fz * (PVY1 + PVY2 * dfz) * LVY * LMUY_prime * zeta2 + SVyg   # (4.E29)

    # Low-speed smoothing
    if ip.reductionSmooth != 1.0
        SHy = SHy * ip.reductionSmooth
        SVy = SVy * ip.reductionSmooth
    end

    alphay = alpha_star + SHy                                   # (4.E20)

    # ── Shape, friction peak, curvature ──────────────────────────────────────
    Cy  = PCY1 * LCY                                            # (4.E21)
    muy = (PDY1 + PDY2 * dfz) *
          (1.0 + PPY3 * dpi + PPY4 * dpi^2) *
          (1.0 - PDY3 * gamma_star^2) * LMUY_star               # (4.E23)

    Dy  = muy * Fz * zeta2                                      # (4.E22)

    signAlphaY = (alphay >= 0.0) ? 1.0 : -1.0
    Ey = (PEY1 + PEY2 * dfz) *
         (1.0 + PEY5 * gamma_star^2 -
         (PEY3 + PEY4 * gamma_star) * signAlphaY) * p.ley       # (4.E24)
    Ey = min(Ey, 1.0)

    signDy = (Dy >= 0.0) ? 1.0 : -1.0
    By     = Kya / (Cy * Dy + epsilony * signDy)                # (4.E26)

    # ── Pure lateral force ────────────────────────────────────────────────────
    Fy0 = Dy * sin(Cy * atan(By * alphay -
          Ey * (By * alphay - atan(By * alphay)))) + SVy         # (4.E19)

    if modes.use_alpha_star && Vcx < 0.0; Fy0 = -Fy0; end

    muy = (Fz == 0.0) ? 0.0 : muy                               # zero-Fz guard

    # Update zeta2 in InternalParams (consumed by Mz0 and Fy combined)
    ip_updated = InternalParams(ip.epsilonx, ip.epsilonk, ip.epsilony,
                                ip.epsilonr, ip.epsilonv, ip.epsilong,
                                ip.reductionSmooth, ip.reductionSharp,
                                ip.reductionLinear, ip.reductionLinear_alpha,
                                zeta2)

    return (Fy0, muy,
            Fy0Intermediates(SHy, SVy, By, Cy, Kya, Kyg0, SVyg),
            ip_updated)
end