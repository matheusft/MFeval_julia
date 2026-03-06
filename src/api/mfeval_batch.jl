# ==============================================================================
# mfeval_batch.jl
#
# Batch entry point: N×6 (or N×7, N×8) input matrix → N×30 output matrix.
# Each row is fully independent so we use Threads.@threads for free
# linear scaling across cores.
#
# Memory strategy
# ---------------
#   • The caller pre-allocates out::Matrix{Float64}(undef, N, 30).
#     Zero allocation per call after the first.
#   • Each thread evaluates one row via the scalar `mfeval` path —
#     the entire solver pipeline is inlined into the loop body by the JIT.
#   • `inst_Kya` (col 29) is computed via finite difference from adjacent
#     rows AFTER the parallel pass, then written back in a single serial pass.
#
# Thread safety
# -------------
#   Row i writes only to out[i, :], so there are no data races.
# ==============================================================================

"""
    mfeval!(out, p, inputs, modes)

In-place batch evaluation.  Writes results into the caller-supplied matrix
`out` (size N×30) to avoid any allocation on the hot path.

# Arguments
- `out    :: Matrix{Float64}` — pre-allocated N×30 output matrix (modified in-place)
- `p      :: TireParams`      — tyre parameter set
- `inputs :: Matrix{Float64}` — N×6, N×7, or N×8 input matrix
  - col 1: Fz [N]
  - col 2: kappa [-]
  - col 3: alpha [rad]
  - col 4: gamma [rad]
  - col 5: phit  [1/m]
  - col 6: Vx    [m/s]
  - col 7: pressure [Pa]  (optional; 0 → use TIR nominal)
  - col 8: omega [rad/s]  (optional; 0 → compute from slip)
- `modes  :: MFModes`         — behaviour flags

# Returns
`out` (the same matrix, for chaining convenience).
"""
function mfeval!(out    ::Matrix{Float64},
                 p      ::TireParams,
                 inputs ::Matrix{Float64},
                 modes  ::MFModes) ::Matrix{Float64}

    N    = size(inputs, 1)
    ncol = size(inputs, 2)

    @assert size(out, 1) == N    "out has $(size(out,1)) rows but inputs has $N"
    @assert size(out, 2) == 30   "out must have 30 columns"
    @assert ncol >= 6             "inputs must have at least 6 columns"

    # ── Parallel evaluation ───────────────────────────────────────────────────
    Threads.@threads for i in 1:N
        pressure = ncol >= 7 ? inputs[i, 7] : 0.0
        omega    = ncol >= 8 ? inputs[i, 8] : 0.0

        inp = MFInputs(inputs[i,1], inputs[i,2], inputs[i,3],
                       inputs[i,4], inputs[i,5], inputs[i,6],
                       pressure, omega)

        r = mfeval(p, inp, modes)

        # Write row in MATLAB column order
        out[i,  1] = r.Fx;       out[i,  2] = r.Fy;       out[i,  3] = r.Fz
        out[i,  4] = r.Mx;       out[i,  5] = r.My;       out[i,  6] = r.Mz
        out[i,  7] = r.kappa;    out[i,  8] = r.alpha;    out[i,  9] = r.gamma
        out[i, 10] = r.phit;     out[i, 11] = r.Vx;       out[i, 12] = r.pressure
        out[i, 13] = r.Re;       out[i, 14] = r.rho;      out[i, 15] = r.two_a
        out[i, 16] = r.t;        out[i, 17] = r.mux;      out[i, 18] = r.muy
        out[i, 19] = r.omega;    out[i, 20] = r.Rl;       out[i, 21] = r.two_b
        out[i, 22] = r.Mzr;     out[i, 23] = r.Cx;       out[i, 24] = r.Cy
        out[i, 25] = r.Cz;      out[i, 26] = r.Kya;      out[i, 27] = r.sigmax
        out[i, 28] = r.sigmay;  out[i, 29] = 0.0;         out[i, 30] = r.Kxk
    end

    # ── inst_Kya (col 29) — finite-difference post-pass (matches MATLAB) ─────
    # MATLAB: diffFY = diff(Fy)./diff(SA); instKya = [diffFY; diffFY(end)]
    if N >= 2
        @inbounds for i in 1:N-1
            d_alpha = out[i+1, 8] - out[i, 8]
            d_Fy    = out[i+1, 2] - out[i, 2] 
            out[i, 29] = abs(d_alpha) < 1e-15 ? 0.0 : d_Fy / d_alpha
        end
        @inbounds out[N, 29] = out[N-1, 29]
    end

    return out
end

"""
    mfeval(p, inputs, modes) → Matrix{Float64}

Allocating variant of the batch evaluation.  Pre-allocates the N×30
output matrix and calls [`mfeval!`](@ref).

Prefer `mfeval!` in performance-critical loops to avoid repeated allocation.

# Arguments
Same as `mfeval!` except no `out` argument.

# Returns
A freshly allocated `N×30 Matrix{Float64}`.
"""
function mfeval(p      ::TireParams,
                inputs ::Matrix{Float64},
                modes  ::MFModes) ::Matrix{Float64}
    N   = size(inputs, 1)
    out = Matrix{Float64}(undef, N, 30)
    return mfeval!(out, p, inputs, modes)
end