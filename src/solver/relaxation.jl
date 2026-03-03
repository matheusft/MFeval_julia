# ==============================================================================
# relaxation.jl
#
# Tyre stiffness and relaxation lengths:
#   calc_relax    — Cx, Cy, sigmax, sigmay
#   calc_inst_Kya — instantaneous cornering stiffness (batch only)
#
# MF5.2 uses the PTX/PTY relaxation-length equations from the MF-Tyre manual.
# MF6.1/6.2 derive relaxation lengths from Kxk/Cx and Kya/Cy (Eqn 19 paper).
#
# Reference: Besselink et al. paper, Eqns (17–19)
#            MF-Tyre equation manual (MF5.2)
# ==============================================================================

"""
    calc_relax(p, pp, iv, Kxk, Kya) → (Cx, Cy, sigmax, sigmay)

Overall longitudinal stiffness (N/m), lateral stiffness (N/m), and
relaxation lengths (m).
"""
@inline @fastmath function calc_relax(p    ::TireParams,
                                       pp   ::PostProInputs,
                                       iv   ::IncrVars,
                                       Kxk  ::Float64,
                                       Kya  ::Float64) ::Tuple{Float64,Float64,Float64,Float64}

    Fz    = pp.Fz
    gamma = pp.gamma
    press = pp.p

    LFZO  = p.lfzo
    Fz0   = p.fnomin
    R0    = p.unloaded_radius
    pi0   = p.nompres
    cx0   = p.longitudinal_stiffness
    cy0   = p.lateral_stiffness
    PCFX1 = p.pcfx1; PCFX2 = p.pcfx2; PCFX3 = p.pcfx3
    PCFY1 = p.pcfy1; PCFY2 = p.pcfy2; PCFY3 = p.pcfy3

    Fz0_prime = LFZO * Fz0
    dfz = (Fz - Fz0_prime) / Fz0_prime
    dpi = (press - pi0) / pi0

    # Overall carcass stiffnesses  Eqns (17) and (18) paper
    Cx = cx0 * (1.0 + PCFX1 * dfz + PCFX2 * dfz^2) * (1.0 + PCFX3 * dpi)
    Cy = cy0 * (1.0 + PCFY1 * dfz + PCFY2 * dfz^2) * (1.0 + PCFY3 * dpi)

    fittyp = p.metadata.fittyp
    sigmax, sigmay = if fittyp == 6 || fittyp == 21
        # MF5.2 — from MF-Tyre equation manual
        PTX1  = p.ptx1; PTX2 = p.ptx2; PTX3 = p.ptx3
        PTY1  = p.pty1; PTY2 = p.pty2
        LSGKP = p.lsgkp; LSGAL = p.lsgal
        PKY3  = p.pky3

        sx = (PTX1 + PTX2 * dfz) * exp(-PTX3 * dfz) * LSGKP * R0 * Fz / Fz0
        sigmayg    = 1.0 - PKY3 * abs(gamma)
        PTYfzn     = PTY2 * Fz0_prime
        sy = PTY1 * sin(2.0 * atan(Fz / PTYfzn)) * sigmayg * R0 * LFZO * LSGAL
        (sx, sy)
    else
        # MF6.1 / MF6.2 — Eqn (19) paper
        (abs(Kxk / Cx), abs(Kya / Cy))
    end

    return (Cx, Cy, sigmax, sigmay)
end

"""
    calc_inst_Kya(alpha, Fy_prev, Fy_curr) → Float64

Instantaneous cornering stiffness ∂Fy/∂α from two adjacent evaluation
points.  Returns `0.0` when called on a single point (batch mode only).

In the single-point API this always returns 0.0; the batch API calls this
after filling the output matrix by differencing adjacent rows.
"""
@inline function calc_inst_Kya(d_alpha ::Float64,
                                d_Fy    ::Float64) ::Float64
    abs(d_alpha) < 1e-15 ? 0.0 : -d_Fy / d_alpha
end