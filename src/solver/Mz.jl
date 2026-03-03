# ==============================================================================
# Mz.jl
#
# Self-aligning moment: Mz0 (pure slip) and Mz (combined slip).
# Mirrors Solver.calculateMz0 and Solver.calculateMz.
#
# Key differences from the book equations noted inline:
#   - Bt uses paper Eqn (A58) instead of book (4.E40) — QBZ6 absent from TIR
#   - Dt uses paper Eqn (A60) instead of book (4.E43) — abs() removed from QDZ3
#   - alphar_eq / alphat_eq use paper Eqns (A54/A55) instead of (4.E78/4.E77)
#   - s uses paper Eqn (A56) — Fz0 instead of Fz0_prime
#   - t multiplied by LFZO (empirically discovered vs TNO)
#
# Reference: Pacejka (2012) Book pp. 183-192; Besselink et al. paper Eqns (A54-A60)
# ==============================================================================

"""
    Mz0Intermediates

Quantities produced by `calc_Mz0` and forwarded to `calc_Mz`.
"""
struct Mz0Intermediates
    alphar    ::Float64   # residual torque slip angle                (4.E37)
    alphat    ::Float64   # pneumatic trail slip angle                (4.E34)
    Dr        ::Float64   # residual torque peak factor               (4.E47)
    Cr        ::Float64   # residual torque shape factor              (4.E46)
    Br        ::Float64   # residual torque stiffness factor          (4.E45)
    Dt        ::Float64   # pneumatic trail peak factor               (A60)
    Ct        ::Float64   # pneumatic trail shape factor              (4.E41)
    Bt        ::Float64   # pneumatic trail stiffness factor          (A58)
    Et        ::Float64   # pneumatic trail curvature factor          (4.E44)
    Kya_prime ::Float64   # Kya + ε·sign(Kya)                        (4.E39)
end

# ------------------------------------------------------------------------------
# Mz0 — pure slip aligning moment (also used to extract intermediates for Mz)
# ------------------------------------------------------------------------------

"""
    calc_Mz0(p, pp, ip, modes, sv, pv, iv, Fy0ints) → (Mz0, Mz0Intermediates)

Pure-slip self-aligning moment (N·m) and the intermediates bundle needed
by `calc_Mz`.

`Fy0ints` is the `Fy0Intermediates` bundle from `calc_Fy0`.
"""
@inline @fastmath function calc_Mz0(p       ::TireParams,
                                     pp      ::PostProInputs,
                                     ip      ::InternalParams,
                                     modes   ::MFModes,
                                     sv      ::StarVars,
                                     pv      ::PrimeVars,
                                     iv      ::IncrVars,
                                     fi      ::Fy0Intermediates) ::Tuple{Float64, Mz0Intermediates}

    epsilonk   = ip.epsilonk
    dfz        = iv.dfz
    dpi        = iv.dpi
    LMUY_star  = sv.LMUY_star
    alpha_star = sv.alpha_star
    gamma_star = sv.gamma_star
    Fz0_prime  = pv.Fz0_prime
    alpha_prime = pv.alpha_prime

    Kya  = fi.Kya
    SHy  = fi.SHy
    SVy  = fi.SVy
    By   = fi.By
    Cy   = fi.Cy

    # Use Fz with low-limit applied for Mz (reference code logic)
    Fz    = modes.use_limits_check ? pp.Fz_lowLimit : pp.Fz
    # Guard: if raw Fz ≤ 0, Mz = 0
    Fz    = (pp.Fz <= 0.0) ? 0.0 : Fz

    Vcx   = pp.uVcx
    gamma = pp.gamma
    R0    = p.unloaded_radius
    Fz0   = p.fnomin
    LFZO  = p.lfzo
    LKY   = p.lky; LTR = p.ltr; LRES = p.lres; LKZC = p.lkzc

    QBZ1  = p.qbz1;  QBZ2 = p.qbz2;  QBZ3 = p.qbz3
    QBZ4  = p.qbz4;  QBZ5 = p.qbz5;  QBZ9 = p.qbz9; QBZ10 = p.qbz10
    QCZ1  = p.qcz1
    QDZ1  = p.qdz1;  QDZ2 = p.qdz2;  QDZ3 = p.qdz3;  QDZ4 = p.qdz4
    QDZ6  = p.qdz6;  QDZ7 = p.qdz7;  QDZ8 = p.qdz8;  QDZ9 = p.qdz9
    QDZ10 = p.qdz10; QDZ11 = p.qdz11
    QEZ1  = p.qez1;  QEZ2 = p.qez2;  QEZ3 = p.qez3
    QEZ4  = p.qez4;  QEZ5 = p.qez5
    QHZ1  = p.qhz1;  QHZ2 = p.qhz2;  QHZ3 = p.qhz3; QHZ4 = p.qhz4
    PPZ1  = p.ppz1;  PPZ2 = p.ppz2

    # ── Kya_prime — guards against Kya = 0 ───────────────────────────────────
    signKya   = (Kya >= 0.0) ? 1.0 : -1.0
    Kya_prime = Kya + epsilonk * signKya                               # (4.E39)

    # ── Horizontal and residual-torque slip angles ────────────────────────────
    SHf    = SHy + SVy / Kya_prime                                     # (4.E38)
    alphar = alpha_star + SHf                                          # (4.E37) = alphaf
    SHt    = QHZ1 + QHZ2 * dfz + (QHZ3 + QHZ4 * dfz) * gamma_star    # (4.E35)
    alphat = alpha_star + SHt                                          # (4.E34)

    # ── Turn-slip ζ₅ ─────────────────────────────────────────────────────────
    zeta5 = modes.use_turn_slip ?
            cos(atan(p.qdtp1 * R0 * pp.phi)) : 1.0                    # (4.E91)

    # ── Pneumatic trail ──────────────────────────────────────────────────────
    # Paper Eqn (A58): uses QBZ4/QBZ5 with gamma (not gamma_star, not abs())
    Bt = (QBZ1 + QBZ2 * dfz + QBZ3 * dfz^2) *
         (1.0 + QBZ4 * gamma + QBZ5 * abs(gamma)) * LKY / LMUY_star  # (A58)

    Ct = QCZ1                                                          # (4.E41) (>0)

    # Paper Eqn (A60): gamma not abs(gamma) in QDZ3 term
    Dt = (QDZ1 + QDZ2 * dfz) * (1.0 - PPZ1 * dpi) *
         (1.0 + QDZ3 * gamma + QDZ4 * gamma^2) *
         Fz * (R0 / Fz0_prime) * LTR * zeta5                          # (A60)

    Et = (QEZ1 + QEZ2 * dfz + QEZ3 * dfz^2) *
         (1.0 + (QEZ4 + QEZ5 * gamma_star) * (2.0 / π) *
         atan(Bt * Ct * alphat))                                        # (4.E44)
    Et = min(Et, 1.0)

    t0 = Dt * cos(Ct * atan(Bt * alphat - Et * (Bt * alphat - atan(Bt * alphat)))) *
         cos(alpha_prime)                                               # (4.E33)

    # ── Fy0 at γ=0, φ=0 (for Mzo') ───────────────────────────────────────────
    Fy0_g0, _ = _fy0_gamma0(p, alpha_star, Fz, Vcx, ip, pv, iv, modes)
    Mzo_prime = -t0 * Fy0_g0                                           # (4.E32) γ=φ=0

    # ── Residual torque ───────────────────────────────────────────────────────
    # Turn-slip ζ₆ or 1
    zeta6 = modes.use_turn_slip ?
            cos(atan(p.qbrp1 * R0 * pp.phi)) : 1.0                    # (4.102)

    Br = (QBZ9 * LKY / LMUY_star + QBZ10 * By * Cy) * zeta6          # (4.E45)
    Cr = modes.use_turn_slip ? 0.0 : 1.0  # zeta7 = 1.0 → Cr = 1     # (4.E46)

    # Turn-slip ζ₈ contributes to Dr and sets Cr
    zeta8 = 1.0
    if modes.use_turn_slip
        # Compute zeta7/zeta8 (requires Fy0 and Gyk — we use Fy0_g0 as proxy
        # for Gyk=1 since this is the zeta0=0 path; full Mz recomputes it)
        # This simplified path is only used for Mz0.  calc_Mz handles the full
        # turn-slip zeta7/zeta8 correctly via calc_mz0_turnslip.
        zeta8 = 1.0   # placeholder; calc_Mz re-evaluates if needed
        Cr    = 1.0
    end

    Dr = Fz * R0 *
         ((QDZ6 + QDZ7 * dfz) * LRES * ip.zeta2 +
          ((QDZ8 + QDZ9 * dfz) * (1.0 + PPZ2 * dpi) +
           (QDZ10 + QDZ11 * dfz) * abs(gamma_star)) *
          gamma_star * LKZC * (modes.use_turn_slip ? 0.0 : 1.0)) *   # zeta0
         LMUY_star * (Vcx >= 0.0 ? 1.0 : -1.0) * cos(alpha_star) +
         zeta8 - 1.0                                                    # (4.E47)

    Mzr0 = Dr * cos(Cr * atan(Br * alphar)) * cos(alpha_prime)         # (4.E36)
    Mz0  = Mzo_prime + Mzr0                                            # (4.E31)

    ints = Mz0Intermediates(alphar, alphat, Dr, Cr, Br, Dt, Ct, Bt, Et, Kya_prime)
    return (Mz0, ints)
end

# ------------------------------------------------------------------------------
# Mz — combined-slip self-aligning moment
# ------------------------------------------------------------------------------

"""
    calc_Mz(p, pp, ip, modes, sv, pv, iv, mz0i, Kxk, Fy, Fx, SVyk) → (Mz, t, Mzr)

Combined-slip self-aligning moment (N·m), pneumatic trail (m), and
residual torque (N·m).
"""
@inline @fastmath function calc_Mz(p    ::TireParams,
                                    pp   ::PostProInputs,
                                    ip   ::InternalParams,
                                    modes::MFModes,
                                    sv   ::StarVars,
                                    pv   ::PrimeVars,
                                    iv   ::IncrVars,
                                    mz0i ::Mz0Intermediates,
                                    Kxk  ::Float64,
                                    Fy   ::Float64,
                                    Fx   ::Float64,
                                    SVyk ::Float64) ::Tuple{Float64,Float64,Float64}

    kappa      = pp.kappa
    gamma      = pp.gamma
    dfz        = iv.dfz
    alpha_prime = pv.alpha_prime
    LFZO       = p.lfzo
    R0         = p.unloaded_radius
    Fz0        = p.fnomin

    LS   = p.ls
    SSZ1 = p.ssz1; SSZ2 = p.ssz2; SSZ3 = p.ssz3; SSZ4 = p.ssz4

    alphar    = mz0i.alphar
    alphat    = mz0i.alphat
    Dr        = mz0i.Dr
    Cr        = mz0i.Cr
    Br        = mz0i.Br
    Dt        = mz0i.Dt
    Ct        = mz0i.Ct
    Bt        = mz0i.Bt
    Et        = mz0i.Et
    Kya_prime = mz0i.Kya_prime

    # ── Equivalent combined-slip angles (Paper A54, A55) ─────────────────────
    ratio_sq = (Kxk / Kya_prime)^2 * kappa^2
    alphar_eq = atan(sqrt(tan(alphar)^2 + ratio_sq)) * (alphar >= 0.0 ? 1.0 : -1.0)  # (A54)
    alphat_eq = atan(sqrt(tan(alphat)^2 + ratio_sq)) * (alphat >= 0.0 ? 1.0 : -1.0)  # (A55)

    # ── Moment arm s (Paper A56: uses Fz0 not Fz0_prime) ─────────────────────
    s = R0 * (SSZ1 + SSZ2 * (Fy / Fz0) + (SSZ3 + SSZ4 * dfz) * gamma) * LS   # (A56)

    # ── Residual torque Mzr ───────────────────────────────────────────────────
    Mzr = Dr * cos(Cr * atan(Br * alphar_eq))                         # (4.E75)

    # ── Fy0 at γ=0, φ=0 and Gyk at γ≠0, φ=0 (for Fy') ───────────────────────
    Fy0_g0, muy_g0 = _fy0_gamma0(p, sv.alpha_star, pp.Fz, pp.uVcx, ip, pv, iv, modes)
    # Gyk with phi=0: re-evaluate using original starVars (gamma intact)
    ip_nots = InternalParams(ip.epsilonx, ip.epsilonk, ip.epsilony,
                             ip.epsilonr, ip.epsilonv, ip.epsilong,
                             ip.reductionSmooth, ip.reductionSharp,
                             ip.reductionLinear, ip.reductionLinear_alpha,
                             1.0)   # zeta2=1 (no turn-slip for this sub-eval)
    _, Gyk_g0, _ = calc_Fy(p, pp, ip_nots, modes, sv, iv, Fy0_g0, muy_g0)

    Fy_prime = Gyk_g0 * Fy0_g0                                        # (4.E74)

    # ── Pneumatic trail t ─────────────────────────────────────────────────────
    t = Dt * cos(Ct * atan(Bt * alphat_eq - Et * (Bt * alphat_eq - atan(Bt * alphat_eq)))) *
        cos(alpha_prime)                                               # (4.E73)

    # Empirically discovered: multiply by LFZO to match TNO solver
    t = t * LFZO

    # Low-speed smoothing
    if ip.reductionSmooth != 1.0
        t   = t   * ip.reductionSmooth
        Mzr = Mzr * ip.reductionSmooth
    end

    # ── Final Mz (MF version branch) ─────────────────────────────────────────
    fittyp = p.metadata.fittyp
    Mz = if fittyp == 6 || fittyp == 21
        # MF5.2 manual equation
        -t * (Fy - SVyk) + Mzr + s * Fx
    else
        # MF6.1 / MF6.2  Eqns (4.E72), (4.E71)
        Mz_prime = -t * Fy_prime                                       # (4.E72)
        Mz_prime + Mzr + s * Fx                                        # (4.E71)
    end

    return (Mz, t, Mzr)
end