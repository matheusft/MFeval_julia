# ==============================================================================
# MFInputs.jl
#
# Single-point input struct for one Magic Formula evaluation.
# All fields are Float64 scalars — the batch layer handles arrays.
# ==============================================================================

"""
    MFInputs(Fz, kappa, alpha, gamma, phit, Vx[, pressure[, omega]])

Single-point input to the Magic Formula solver. All quantities in SI units.

Fields
------
- `Fz`       : Normal load on the tyre                   [N]
- `kappa`    : Longitudinal slip (−1 = locked wheel)     [−]
- `alpha`    : Side-slip angle                            [rad]
- `gamma`    : Inclination (camber) angle                 [rad]
- `phit`     : Turn slip                                  [1/m]
- `Vx`       : Forward (longitudinal) velocity            [m/s]
- `pressure` : Inflation pressure (0 → use TIR INFLPRES)  [Pa]
- `omega`    : Wheel rotational speed (0 → estimated)    [rad/s]
"""
struct MFInputs
    Fz       ::Float64
    kappa    ::Float64
    alpha    ::Float64
    gamma    ::Float64
    phit     ::Float64
    Vx       ::Float64
    pressure ::Float64   # 0.0 → fall back to TIR INFLPRES
    omega    ::Float64   # 0.0 → estimated internally from Re
end

# Convenience constructors for 6-column and 7-column input matrices
MFInputs(Fz, kappa, alpha, gamma, phit, Vx) =
    MFInputs(Fz, kappa, alpha, gamma, phit, Vx, 0.0, 0.0)

MFInputs(Fz, kappa, alpha, gamma, phit, Vx, pressure) =
    MFInputs(Fz, kappa, alpha, gamma, phit, Vx, pressure, 0.0)