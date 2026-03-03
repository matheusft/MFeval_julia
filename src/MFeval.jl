module MFeval

# ── Phase 1: types ─────────────────────────────────────────────────────────────
include("types/TireParams.jl")
include("types/MFInputs.jl")
include("types/MFModes.jl")
include("types/MFOutputs.jl")

# ── Phase 1: I/O ───────────────────────────────────────────────────────────────
include("io/read_tir.jl")

# ── Phase 2: solver (uncomment as each module is implemented) ──────────────────
# include("solver/parse_inputs.jl")
# include("solver/basic_vars.jl")
# include("solver/Fx0.jl")
# include("solver/Fy0.jl")
# include("solver/combined_slip.jl")
# include("solver/Mx.jl")
# include("solver/My.jl")
# include("solver/Mz.jl")
# include("solver/turn_slip.jl")
# include("solver/geometry.jl")
# include("solver/relaxation.jl")

# ── Phase 3: public API (uncomment when solver is complete) ────────────────────
# include("api/mfeval.jl")
# include("api/mfeval_batch.jl")

# ── Exports ────────────────────────────────────────────────────────────────────
export TireParams, TireMetadata
export MFInputs
export MFModes
export MFOutputs, to_vector, from_matrix_row
export read_tir
export version, is_mf52, is_mf61, is_mf62

end # module MFeval