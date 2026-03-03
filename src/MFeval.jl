module MFeval

# ── Phase 1: types ─────────────────────────────────────────────────────────────
include("types/TireParams.jl")
include("types/MFInputs.jl")
include("types/MFModes.jl")
include("types/MFOutputs.jl")

# ── Phase 1: I/O ───────────────────────────────────────────────────────────────
include("io/read_tir.jl")

# ── Phase 2: solver ────────────────────────────────────────────────────────────
# Dependency order matters: each file may reference types defined above it.
include("solver/parse_inputs.jl")    # PostProInputs, InternalParams
include("solver/basic_vars.jl")      # StarVars, PrimeVars, IncrVars, SlipVelocities
include("solver/turn_slip.jl")       # TurnSlipFactors, calc_zeta*
include("solver/Fx0.jl")             # calc_Fx0
include("solver/Fy0.jl")             # calc_Fy0, Fy0Intermediates
include("solver/combined_slip.jl")   # calc_Fx, calc_Fy
include("solver/Mx.jl")              # calc_Mx
include("solver/My.jl")              # calc_My
include("solver/Mz.jl")              # calc_Mz0, calc_Mz, Mz0Intermediates
include("solver/geometry.jl")        # calc_Re, calc_rho_Rl, calc_contact_patch
include("solver/relaxation.jl")      # calc_relax, calc_inst_Kya

# ── Phase 3: public API (uncomment when complete) ──────────────────────────────
# include("api/mfeval.jl")
# include("api/mfeval_batch.jl")

# ── Exports ────────────────────────────────────────────────────────────────────
export TireParams, TireMetadata
export MFInputs
export MFModes
export MFOutputs, to_vector, from_matrix_row
export read_tir
export version, is_mf52, is_mf61, is_mf62

# Phase 2 solver internals (exported for testing and power-user access)
export PostProInputs, InternalParams
export StarVars, PrimeVars, IncrVars, SlipVelocities
export Fy0Intermediates, Mz0Intermediates
export TurnSlipFactors
export calc_zeta1, calc_zeta2_zeta3, calc_mz0_turnslip
export parse_inputs
export calc_basic_vars
export calc_Fx0, calc_Fy0, calc_Fx, calc_Fy
export calc_Mx, calc_My
export calc_Mz0, calc_Mz
export calc_Re, calc_rho_Rl, calc_contact_patch
export calc_relax, calc_inst_Kya

end # module MFeval