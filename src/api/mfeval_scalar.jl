# ==============================================================================
# mfeval_scalar.jl
#
# Public scalar entry point: one MFInputs → one MFOutputs.
# Calls the full Phase 2 pipeline in dependency order, then packs all
# 30 output columns into an MFOutputs struct.
#
# Design constraints (Phase 5 targets preserved):
#   • @inline on the hot path — lets the batch loop inline the whole pipeline
#   • @fastmath throughout — safe for tyre model arithmetic
#   • Zero heap allocations — all intermediates are stack-allocated structs
#   • No runtime FITTYP checks — dispatch via TireParams{V}
# ==============================================================================

"""
    mfeval(p, inp, modes) → MFOutputs

Evaluate the Pacejka Magic Formula for a single operating point.

# Arguments
- `p      :: TireParams`  — tyre parameter set from `read_tir`
- `inp    :: MFInputs`    — scalar input (Fz, κ, α, γ, φ, Vx [, pressure [, ω]])
- `modes  :: MFModes`     — behaviour flags decoded from a 3-digit `useMode` integer

# Returns
An [`MFOutputs`](@ref) struct with all 30 output quantities.

# Example
```julia
p    = read_tir("MagicFormula61_Parameters.tir")
inp  = MFInputs(4000.0, 0.0, 0.05, 0.0, 0.0, 16.7)
m    = MFModes(111)
out  = mfeval(p, inp, m)
println(out.Fy)   # lateral force [N]
```
"""
@inline @fastmath function mfeval(p     ::TireParams,
                                   inp   ::MFInputs,
                                   modes ::MFModes) ::MFOutputs

    # Check if useMode is valid (matches MATLAB behavior)
    if !modes.is_valid
        # Return NaN for all outputs, matching MATLAB behavior
        return MFOutputs(
            NaN, NaN, NaN, NaN, NaN, NaN,                    # Forces and moments
            NaN, NaN, NaN, NaN, NaN, NaN,                    # Input echoes
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, # Geometry
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN           # Stiffnesses
        )
    end

    # ── Stage 1: pre-processing ───────────────────────────────────────────────
    pp, ip = parse_inputs(p, inp, modes)

    # ── Stage 2: shared variables ─────────────────────────────────────────────
    sv, pv, iv, _ = calc_basic_vars(p, pp, ip, modes)

    # ── Stage 3: pure-slip forces ─────────────────────────────────────────────
    Fx0, mux, Kxk          = calc_Fx0(p, pp, ip, modes, sv, pv, iv)
    Fy0, muy, fi, ip2      = calc_Fy0(p, pp, ip, modes, sv, pv, iv)

    # ── Stage 4: aligning torque intermediates ────────────────────────────────
    Mz0, mz0i = calc_Mz0(p, pp, ip2, modes, sv, pv, iv, fi)

    # ── Stage 5: combined-slip forces ─────────────────────────────────────────
    Fx              = calc_Fx(p, pp, modes, sv, iv, Fx0)
    Fy, _Gyk, SVyk  = calc_Fy(p, pp, ip2, modes, sv, iv, Fy0, muy)

    # ── Stage 6: moments ──────────────────────────────────────────────────────
    Mx             = calc_Mx(p, pp, iv, Fy)
    My             = calc_My(p, pp, Fx)
    Mz, t, Mzr    = calc_Mz(p, pp, ip2, modes, sv, pv, iv, mz0i, Kxk, Fy, Fx, SVyk)

    # ── Stage 7: geometry ─────────────────────────────────────────────────────
    dpi            = iv.dpi
    Re, Romega, ω  = calc_Re(p, pp, dpi)
    rho, Rl, Cz    = calc_rho_Rl(p, pp, dpi, ω, Romega, Fx, Fy)
    a, b, _NCz     = calc_contact_patch(p, pp, dpi)

    # ── Stage 8: relaxation lengths ───────────────────────────────────────────
    Cx, Cy, σx, σy = calc_relax(p, pp, iv, Kxk, fi.Kya)

    # ── Pack outputs (column order matches MATLAB mfeval exactly) ─────────────
    return MFOutputs(
        Fx, Fy, pp.uFz, Mx, My, Mz,        # cols  1-6
        pp.ukappaLow, pp.ualpha, pp.ugamma, # cols  7-9
        pp.uphit, pp.uVcx, pp.p,            # cols 10-12
        Re, rho, 2.0*a, t,                  # cols 13-16
        mux, muy, ω, Rl, 2.0*b, Mzr,       # cols 17-22
        Cx, Cy, Cz, fi.Kya, σx, σy,        # cols 23-28
        0.0,                                # col  29: inst_Kya — N/A for scalar call
        Kxk,                                # col  30
    )
end