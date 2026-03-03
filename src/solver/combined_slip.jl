# ==============================================================================
# combined_slip.jl
#
# Combined-slip weighting functions:
#   Fx = Gxa · Fx0            Eqn (4.E50)
#   Fy = Gyk · Fy0 + SVyk     Eqn (4.E58)
#
# Mirrors Solver.calculateFx and Solver.calculateFy.
# Reference: Pacejka (2012) Book, pp. 187-192 (4.E50–4.E67)
# ==============================================================================

"""
    calc_Fx(p, pp, modes, sv, iv, Fx0) → Fx

Combined longitudinal force (N).
"""
@inline @fastmath function calc_Fx(p     ::TireParams,
                                    pp    ::PostProInputs,
                                    modes ::MFModes,
                                    sv    ::StarVars,
                                    iv    ::IncrVars,
                                    Fx0   ::Float64) ::Float64

    kappa      = pp.kappa
    alpha_star = sv.alpha_star
    gamma_star = sv.gamma_star
    dfz        = iv.dfz

    LXAL = p.lxal
    RBX1 = p.rbx1;  RBX2 = p.rbx2;  RBX3 = p.rbx3
    RCX1 = p.rcx1
    REX1 = p.rex1;  REX2 = p.rex2
    RHX1 = p.rhx1

    Cxa = RCX1                                                          # (4.E55)
    Exa = min(REX1 + REX2 * dfz, 1.0)                                  # (4.E56) ≤1

    SHxa   = RHX1                                                       # (4.E57)
    Bxa    = (RBX1 + RBX3 * gamma_star^2) *
             cos(atan(RBX2 * kappa)) * LXAL                            # (4.E54)
    alphas = alpha_star + SHxa                                          # (4.E53)

    Gxa0 = cos(Cxa * atan(Bxa * SHxa  - Exa * (Bxa * SHxa  - atan(Bxa * SHxa))))   # (4.E52)
    Gxa  = cos(Cxa * atan(Bxa * alphas - Exa * (Bxa * alphas - atan(Bxa * alphas)))) / Gxa0  # (4.E51)

    return Gxa * Fx0                                                    # (4.E50)
end

# ------------------------------------------------------------------------------

"""
    calc_Fy(p, pp, ip, modes, sv, iv, Fy0, muy) → (Fy, Gyk, SVyk)

Combined lateral force (N), weighting factor (-), and kappa-induced
lateral force shift SVyk (N).
"""
@inline @fastmath function calc_Fy(p     ::TireParams,
                                    pp    ::PostProInputs,
                                    ip    ::InternalParams,
                                    modes ::MFModes,
                                    sv    ::StarVars,
                                    iv    ::IncrVars,
                                    Fy0   ::Float64,
                                    muy   ::Float64) ::Tuple{Float64,Float64,Float64}

    Fz         = pp.Fz
    kappa      = pp.kappa
    alpha_star = sv.alpha_star
    gamma_star = sv.gamma_star
    dfz        = iv.dfz
    zeta2      = ip.zeta2

    LYKA  = p.lyka
    LVYKA = p.lvyka
    RBY1  = p.rby1;  RBY2 = p.rby2;  RBY3 = p.rby3;  RBY4 = p.rby4
    RCY1  = p.rcy1
    REY1  = p.rey1;  REY2 = p.rey2
    RHY1  = p.rhy1;  RHY2 = p.rhy2
    RVY1  = p.rvy1;  RVY2 = p.rvy2;  RVY3 = p.rvy3
    RVY4  = p.rvy4;  RVY5 = p.rvy5;  RVY6 = p.rvy6

    # kappa-induced lateral force (Eqn 4.E67)
    DVyk = muy * Fz * (RVY1 + RVY2 * dfz + RVY3 * gamma_star) *
           cos(atan(RVY4 * alpha_star)) * zeta2                         # (4.E67)
    SVyk = DVyk * sin(RVY5 * atan(RVY6 * kappa)) * LVYKA                # (4.E66)

    SHyk = RHY1 + RHY2 * dfz                                            # (4.E65)
    Eyk  = min(REY1 + REY2 * dfz, 1.0)                                  # (4.E64) ≤1

    Cyk    = RCY1                                                        # (4.E63)
    Byk    = (RBY1 + RBY4 * gamma_star^2) *
             cos(atan(RBY2 * (alpha_star - RBY3))) * LYKA               # (4.E62)
    kappas = kappa + SHyk                                               # (4.E61)

    Gyk0 = cos(Cyk * atan(Byk * SHyk  - Eyk * (Byk * SHyk  - atan(Byk * SHyk))))    # (4.E60)
    Gyk  = cos(Cyk * atan(Byk * kappas - Eyk * (Byk * kappas - atan(Byk * kappas)))) / Gyk0  # (4.E59)

    # Low-speed smoothing of SVyk
    if ip.reductionSmooth != 1.0
        SVyk = SVyk * ip.reductionSmooth
    end

    Fy = Gyk * Fy0 + SVyk                                              # (4.E58)

    return (Fy, Gyk, SVyk)
end