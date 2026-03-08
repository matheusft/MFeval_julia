# ==============================================================================
# MFModes.jl
#
# Decoded operating-mode flags for the Magic Formula solver.
# The 3-digit useMode integer encodes three independent boolean switches:
#
#   Hundreds digit  1 → include limit checks      2 → skip limit checks
#   Tens digit      1 → standard alpha            2 → alpha_star
#   Units digit     1 → no turn slip              2 → include turn slip
#
# Decoding once at construction time keeps the hot solver path free of
# integer arithmetic and conditional branches.
# ==============================================================================

"""
    MFModes(use_limits_check, use_alpha_star, use_turn_slip, is_valid)

Operating-mode flags for the Magic Formula solver, decoded from the
3-digit `useMode` integer used by the original MATLAB mfeval.

Use `MFModes(useMode::Integer)` to construct from an integer.
"""
struct MFModes
    use_limits_check ::Bool   # hundreds digit: 1=true, 2=false
    use_alpha_star   ::Bool   # tens digit:     1=false, 2=true
    use_turn_slip    ::Bool   # units digit:    1=false, 2=true
    is_valid         ::Bool   # true if useMode is supported, false → return NaN (matches MATLAB)
end

"""
    MFModes(useMode::Integer) → MFModes

Decode a 3-digit `useMode` integer into an `MFModes` struct.

# Examples
```julia
MFModes(111)  # limits=true,  alpha_star=false, turn_slip=false
MFModes(121)  # limits=true,  alpha_star=true,  turn_slip=false
MFModes(112)  # limits=true,  alpha_star=false, turn_slip=true
MFModes(222)  # limits=false, alpha_star=true,  turn_slip=true
```
"""
function MFModes(useMode::Integer)
    hundreds, rem1 = divrem(useMode, 100)
    tens,     units = divrem(rem1,   10)

    # Check if this is a supported useMode (matches MATLAB behavior)
    is_valid = (hundreds in [1, 2]) && (tens in [1, 2]) && (units in [1, 2])
    
    # For invalid modes, use default values but mark as invalid
    # This matches MATLAB's behavior of accepting the call but returning NaN
    use_limits_check = (hundreds == 1) || !is_valid  # default to true for invalid modes
    use_alpha_star   = (tens == 2) && is_valid       # default to false for invalid modes 
    use_turn_slip    = (units == 2) && is_valid      # default to false for invalid modes

    MFModes(use_limits_check, use_alpha_star, use_turn_slip, is_valid)
end