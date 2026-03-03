# ==============================================================================
# Fx0.jl
#
# Pure longitudinal force Fx0 and longitudinal slip stiffness Kxk.
# Mirrors Solver.calculateFx0 from the reference MATLAB code.
#
# Reference equations: Pacejka (2012) Book, pp. 177-181 (4.E9–4.E19)
# ==============================================================================

"""
    calc_Fx0(p, pp, ip, modes, sv, pv, iv) → (Fx0, mux, Kxk)

Pure longitudinal force (N), longitudinal friction coefficient (-), and
longitudinal slip stiffness (N/-).

Handles turn-slip (ζ₁), low-speed smoothing, and backward-speed correction.
"""
@inline @fastmath function calc_Fx0(p     ::TireParams,
                                     pp    ::PostProInputs,
                                     ip    ::InternalParams,
                                     modes ::MFModes,
                                     sv    ::StarVars,
                                     pv    ::PrimeVars,
                                     iv    ::IncrVars) ::Tuple{Float64,Float64,Float64}

    # ── Scaling coefficients ─────────────────────────────────────────────────
    LCX  = p.lcx
    LEX  = p.lex
    LKX  = p.lkx
    LHX  = p.lhx
    LVX  = p.lvx

    # ── Longitudinal coefficients ─────────────────────────────────────────────
    PCX1 = p.pcx1
    PDX1 = p.pdx1
    PDX2 = p.pdx2
    PDX3 = p.pdx3
    PEX1 = p.pex1
    PEX2 = p.pex2
    PEX3 = p.pex3
    PEX4 = p.pex4
    PKX1 = p.pkx1
    PKX2 = p.pkx2
    PKX3 = p.pkx3
    PHX1 = p.phx1
    PHX2 = p.phx2
    PVX1 = p.pvx1
    PVX2 = p.pvx2
    PPX1 = p.ppx1
    PPX2 = p.ppx2
    PPX3 = p.ppx3
    PPX4 = p.ppx4

    Fz         = pp.Fz
    kappa      = pp.kappa
    gamma      = pp.gamma
    Vcx        = pp.uVcx
    dfz        = iv.dfz
    dpi        = iv.dpi
    LMUX_star  = sv.LMUX_star
    LMUX_prime = pv.LMUX_prime
    epsilonx   = ip.epsilonx

    # ── Turn-slip ζ₁ ─────────────────────────────────────────────────────────
    zeta1 = modes.use_turn_slip ? calc_zeta1(p, pp, dfz) : 1.0

    # ── Shape and peak ───────────────────────────────────────────────────────
    Cx  = PCX1 * LCX                                                      # (4.E11)
    mux = (PDX1 + PDX2 * dfz) *
          (1.0 + PPX3 * dpi + PPX4 * dpi^2) *
          (1.0 - PDX3 * gamma^2) * LMUX_star                             # (4.E13)
    mux = (Fz == 0.0) ? 0.0 : mux                                        # zero-Fz guard
    Dx  = mux * Fz * zeta1                                                # (4.E12)

    # ── Slip stiffness ───────────────────────────────────────────────────────
    Kxk = Fz * (PKX1 + PKX2 * dfz) * exp(PKX3 * dfz) *
          (1.0 + PPX1 * dpi + PPX2 * dpi^2) * LKX                       # (4.E15)

    # ── Stiffness factor Bx ───────────────────────────────────────────────────
    signDx = (Dx >= 0.0) ? 1.0 : -1.0      # sign(0) → 1 to avoid Kxk/0
    Bx  = Kxk / (Cx * Dx + epsilonx * signDx)                           # (4.E16)

    # ── Shifts ───────────────────────────────────────────────────────────────
    SHx = (PHX1 + PHX2 * dfz) * LHX                                     # (4.E17)
    SVx = Fz * (PVX1 + PVX2 * dfz) * LVX * LMUX_prime * zeta1          # (4.E18)

    # Low-speed model: smooth reduction of shifts
    if ip.reductionSmooth != 1.0
        SHx = SHx * ip.reductionSmooth
        SVx = SVx * ip.reductionSmooth
    end

    kappax = kappa + SHx                                                  # (4.E10)

    # ── Curvature factor Ex ───────────────────────────────────────────────────
    Ex = (PEX1 + PEX2 * dfz + PEX3 * dfz^2) *
         (1.0 - PEX4 * (kappax >= 0.0 ? 1.0 : -1.0)) * LEX             # (4.E14)
    Ex = min(Ex, 1.0)   # Ex ≤ 1 (limit check)

    # ── Pure longitudinal force ───────────────────────────────────────────────
    Fx0 = Dx * sin(Cx * atan(Bx * kappax - Ex * (Bx * kappax - atan(Bx * kappax)))) + SVx  # (4.E9)

    # Backward speed: flip sign
    if Vcx < 0.0; Fx0 = -Fx0; end

    return (Fx0, mux, Kxk)
end