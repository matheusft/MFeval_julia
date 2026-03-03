# ==============================================================================
# MFOutputs.jl
#
# Output struct for a single Magic Formula evaluation.
# 30 fields — column order matches the original MATLAB mfeval output matrix
# exactly, so index-based interop and regression tests are unambiguous.
# ==============================================================================

"""
    MFOutputs

All 30 outputs produced by a single Magic Formula evaluation.
Column ordering matches the MATLAB mfeval output matrix exactly.

Forces & moments (cols 1-6)
----------------------------
| Col | Name  | Description                    | Unit  |
|-----|-------|--------------------------------|-------|
|  1  | Fx    | Longitudinal force             | N     |
|  2  | Fy    | Lateral force                  | N     |
|  3  | Fz    | Normal force (pass-through)    | N     |
|  4  | Mx    | Overturning moment             | N·m   |
|  5  | My    | Rolling resistance moment      | N·m   |
|  6  | Mz    | Self-aligning moment           | N·m   |

Slip quantities — post-processed / clamped (cols 7-12)
-------------------------------------------------------
|  7  | kappa    | Longitudinal slip           | −     |
|  8  | alpha    | Side-slip angle             | rad   |
|  9  | gamma    | Inclination angle           | rad   |
| 10  | phit     | Turn slip                   | 1/m   |
| 11  | Vx       | Longitudinal velocity       | m/s   |
| 12  | pressure | Inflation pressure          | Pa    |

Tyre geometry & state (cols 13-22)
------------------------------------
| 13  | Re      | Effective rolling radius     | m     |
| 14  | rho     | Tyre deflection              | m     |
| 15  | two_a   | Contact patch length (2a)    | m     |
| 16  | t       | Pneumatic trail              | m     |
| 17  | mux     | Longitudinal friction coeff. | −     |
| 18  | muy     | Lateral friction coeff.      | −     |
| 19  | omega   | Rotational speed             | rad/s |
| 20  | Rl      | Loaded radius                | m     |
| 21  | two_b   | Contact patch width (2b)     | m     |
| 22  | Mzr     | Residual torque              | N·m   |

Stiffness & relaxation (cols 23-30)
--------------------------------------
| 23  | Cx       | Longitudinal tyre stiffness  | N/m   |
| 24  | Cy       | Lateral tyre stiffness       | N/m   |
| 25  | Cz       | Vertical tyre stiffness      | N/m   |
| 26  | Kya      | Cornering stiffness          | N/rad |
| 27  | sigmax   | Long. relaxation length      | m     |
| 28  | sigmay   | Lat. relaxation length       | m     |
| 29  | inst_Kya | Instantaneous cornering stiff| N/rad |
| 30  | Kxk      | Longitudinal slip stiffness  | N/−   |
"""
struct MFOutputs
    # Forces & moments
    Fx       ::Float64
    Fy       ::Float64
    Fz       ::Float64
    Mx       ::Float64
    My       ::Float64
    Mz       ::Float64
    # Slip (post-processed)
    kappa    ::Float64
    alpha    ::Float64
    gamma    ::Float64
    phit     ::Float64
    Vx       ::Float64
    pressure ::Float64
    # Tyre geometry
    Re       ::Float64
    rho      ::Float64
    two_a    ::Float64
    t        ::Float64
    mux      ::Float64
    muy      ::Float64
    omega    ::Float64
    Rl       ::Float64
    two_b    ::Float64
    Mzr      ::Float64
    # Stiffness & relaxation
    Cx       ::Float64
    Cy       ::Float64
    Cz       ::Float64
    Kya      ::Float64
    sigmax   ::Float64
    sigmay   ::Float64
    inst_Kya ::Float64
    Kxk      ::Float64
end

"""
    to_vector(out::MFOutputs) → Vector{Float64}

Convert an `MFOutputs` to a plain 30-element `Vector{Float64}` in the
original MATLAB column order. Useful for regression tests and interop.
"""
function to_vector(o::MFOutputs) ::Vector{Float64}
    [o.Fx, o.Fy, o.Fz, o.Mx, o.My, o.Mz,
     o.kappa, o.alpha, o.gamma, o.phit, o.Vx, o.pressure,
     o.Re, o.rho, o.two_a, o.t, o.mux, o.muy, o.omega, o.Rl, o.two_b, o.Mzr,
     o.Cx, o.Cy, o.Cz, o.Kya, o.sigmax, o.sigmay, o.inst_Kya, o.Kxk]
end

"""
    from_matrix_row(m::Matrix{Float64}, i::Integer) → MFOutputs

Reconstruct an `MFOutputs` from row `i` of a pre-allocated output matrix.
"""
function from_matrix_row(m::Matrix{Float64}, i::Integer) ::MFOutputs
    MFOutputs(
        m[i, 1],  m[i, 2],  m[i, 3],  m[i, 4],  m[i, 5],  m[i, 6],
        m[i, 7],  m[i, 8],  m[i, 9],  m[i,10],  m[i,11],  m[i,12],
        m[i,13],  m[i,14],  m[i,15],  m[i,16],  m[i,17],  m[i,18],
        m[i,19],  m[i,20],  m[i,21],  m[i,22],
        m[i,23],  m[i,24],  m[i,25],  m[i,26],
        m[i,27],  m[i,28],  m[i,29],  m[i,30],
    )
end