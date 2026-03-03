# ==============================================================================
# Mx.jl
#
# Overturning moment Mx.
# Mirrors Solver.calculateMx from the reference MATLAB code.
#
# Two formulations coexist in the codebase (book 4.E69 vs TNO manual).
# The TNO manual definition is used to match the reference solver output.
#
# Reference: Besselink et al. draft paper, Eqn (49)
#            Pacejka (2012) Book Eqn (4.E69)
# ==============================================================================

"""
    calc_Mx(p, pp, iv, Fy) → Mx

Overturning moment (N·m).

Applies a low-Fz reduction factor empirically present in the reference code:
  Fz_eff = Fz · (Fz/FZMIN)²   when Fz < FZMIN.
"""
@inline @fastmath function calc_Mx(p  ::TireParams,
                                    pp ::PostProInputs,
                                    iv ::IncrVars,
                                    Fy ::Float64) ::Float64

    Fz    = pp.Fz
    gamma = pp.gamma
    dpi   = iv.dpi

    # Low-Fz empirical correction (matches reference MATLAB)
    Fz_eff = (Fz < p.fzmin) ? Fz * (Fz / p.fzmin)^2 : Fz

    R0   = p.unloaded_radius
    Fz0  = p.fnomin
    LVMX = p.lvmx
    LMX  = p.lmx

    QSX1  = p.qsx1;  QSX2  = p.qsx2;  QSX3  = p.qsx3
    QSX4  = p.qsx4;  QSX5  = p.qsx5;  QSX6  = p.qsx6
    QSX7  = p.qsx7;  QSX8  = p.qsx8;  QSX9  = p.qsx9
    QSX10 = p.qsx10; QSX11 = p.qsx11; QSX12 = p.qsx12
    QSX13 = p.qsx13; QSX14 = p.qsx14
    PPMX1 = p.ppmx1

    # TNO Equation Manual definition (matches reference solver)
    # This is also consistent with the draft Besselink paper, Eqn (49)
    Mx = R0 * Fz_eff * LMX *
         (QSX1 * LVMX -
          QSX2 * gamma * (1.0 + PPMX1 * dpi) +
          QSX3 * (Fy / Fz0) +
          QSX4 * cos(QSX5 * atan((QSX6 * Fz_eff / Fz0)^2)) *
                 sin(QSX7 * gamma + QSX8 * atan(QSX9 * Fy / Fz0)) +
          QSX10 * atan(QSX11 * Fz_eff / Fz0) * gamma) +
         R0 * LMX *
         (Fy * (QSX13 + QSX14 * abs(gamma)) -
          Fz_eff * QSX12 * gamma * abs(gamma))

    return Mx
end