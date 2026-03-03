# ==============================================================================
# read_tir.jl
#
# Parser for Tyre Property Files (.tir).
#
# The TIR format is a plain-text key=value file with:
#   - Section headers in square brackets: [SECTION_NAME]
#   - Parameters as:  KEY = value   $optional inline comment
#   - Comment lines starting with '$' or '!'
#   - One section, [SHAPE], that is obsolete and must be skipped entirely
#
# The parser builds a Dict{String,Any}, then passes it to the TireParams
# keyword constructor so that every missing key gets its default value.
# No runtime FITTYP branching is needed here — TireParams handles that.
# ==============================================================================

"""
    read_tir(path::AbstractString) → TireParams

Parse a Magic Formula tyre property file (`.tir`) and return a fully
constructed, version-typed `TireParams` struct.

Throws an informative error if the file cannot be opened or if `FITTYP`
is missing or unsupported.
"""
function read_tir(path::AbstractString) ::TireParams
    isfile(path) || error("read_tir: file not found: \"$path\"")

    raw = _parse_tir_file(path)

    _build_tireparams(raw)
end

# ------------------------------------------------------------------------------
# Internal: parse the file into a flat Dict{String, Union{Float64, String}}
# ------------------------------------------------------------------------------
function _parse_tir_file(path::AbstractString) ::Dict{String, Any}
    dict         = Dict{String, Any}()
    skip_section = false
    mass_count   = 0          # MASS appears twice; first = units string

    open(path, "r") do io
        for line in eachline(io)
            line = strip(line)
            isempty(line) && continue

            # Section header
            if startswith(line, '[')
                skip_section = (line == "[SHAPE]")
                continue
            end

            # Skip comment lines and suppressed section contents
            (startswith(line, '$') || startswith(line, '!') || skip_section) && continue

            # Find the '=' separator
            eq_idx = findfirst('=', line)
            isnothing(eq_idx) && continue

            key = String(strip(line[1:eq_idx-1]))
            isempty(key) && continue

            # Strip inline comment (everything from '$' onward)
            value_raw = line[eq_idx+1:end]
            dollar_idx = findfirst('$', value_raw)
            value_str  = String(strip(isnothing(dollar_idx) ? value_raw : value_raw[1:dollar_idx-1]))
            isempty(value_str) && continue

            # Special-case: MASS appears as a units string first, then a
            # numeric inertia value. We rename the numeric one to MASS1
            # to match the readTIR.m convention.
            if key == "MASS"
                mass_count += 1
                if mass_count == 1
                    # First occurrence → units string (e.g. 'kg')
                    dict["MASS"] = _strip_quotes(value_str)
                else
                    # Second occurrence → inertia value [kg]
                    dict["MASS1"] = _try_parse_float(value_str, key)
                end
                continue
            end

            # Determine whether the value is a number or a string
            dict[key] = _try_parse_float_or_string(value_str, key)
        end
    end

    return dict
end

# ------------------------------------------------------------------------------
# Internal: construct TireParams from the parsed Dict
# ------------------------------------------------------------------------------
function _build_tireparams(d::Dict{String, Any}) ::TireParams

    # FITTYP is mandatory — everything else has a default
    haskey(d, "FITTYP") || error("read_tir: FITTYP not found in tyre property file")
    fittyp = Int(d["FITTYP"])

    meta = TireMetadata(
        tyreside             = _get_str(d, "TYRESIDE",             "LEFT"),
        file_type            = _get_str(d, "FILE_TYPE",            "default"),
        file_version         = _get_f64(d, "FILE_VERSION",          3.0),
        file_format          = _get_str(d, "FILE_FORMAT",          "default"),
        property_file_format = _get_str(d, "PROPERTY_FILE_FORMAT", "default"),
        function_name        = _get_str(d, "FUNCTION_NAME",        "default"),
        units_length         = _get_str(d, "LENGTH",               "meter"),
        units_force          = _get_str(d, "FORCE",                "newton"),
        units_angle          = _get_str(d, "ANGLE",                "radians"),
        fittyp               = fittyp,
        use_mode             = Int(_get_f64(d, "USE_MODE",          124.0)),
        n_tire_states        = Int(_get_f64(d, "N_TIRE_STATES",       4.0)),
        user_sub_id          = Int(_get_f64(d, "USER_SUB_ID",        815.0)),
    )

    return TireParams(
        metadata                = meta,
        # [MODEL]
        longvl                  = _get_f64(d, "LONGVL",              16.7),
        vxlow                   = _get_f64(d, "VXLOW",                1.0),
        road_increment          = _get_f64(d, "ROAD_INCREMENT",       0.01),
        road_direction          = _get_f64(d, "ROAD_DIRECTION",       1.0),
        # [DIMENSION]
        unloaded_radius         = _get_f64(d, "UNLOADED_RADIUS",     0.3135),
        width                   = _get_f64(d, "WIDTH",               0.205),
        aspect_ratio            = _get_f64(d, "ASPECT_RATIO",        0.6),
        rim_radius              = _get_f64(d, "RIM_RADIUS",          0.1905),
        rim_width               = _get_f64(d, "RIM_WIDTH",           0.152),
        # [OPERATING_CONDITIONS]
        inflpres                = _get_f64(d, "INFLPRES",            220000.0),
        nompres                 = _get_f64(d, "NOMPRES",             220000.0),
        # [INERTIA]
        mass                    = _get_f64(d, "MASS1",               9.3),
        ixx                     = _get_f64(d, "IXX",                 0.391),
        iyy                     = _get_f64(d, "IYY",                 0.736),
        belt_mass               = _get_f64(d, "BELT_MASS",           7.0),
        belt_ixx                = _get_f64(d, "BELT_IXX",            0.34),
        belt_iyy                = _get_f64(d, "BELT_IYY",            0.6),
        gravity                 = _get_f64(d, "GRAVITY",             -9.81),
        # [VERTICAL]
        fnomin                  = _get_f64(d, "FNOMIN",              4000.0),
        vertical_stiffness      = _get_f64(d, "VERTICAL_STIFFNESS",  200000.0),
        vertical_damping        = _get_f64(d, "VERTICAL_DAMPING",    50.0),
        mc_contour_a            = _get_f64(d, "MC_CONTOUR_A",        0.5),
        mc_contour_b            = _get_f64(d, "MC_CONTOUR_B",        0.5),
        breff                   = _get_f64(d, "BREFF",               8.4),
        dreff                   = _get_f64(d, "DREFF",               0.27),
        freff                   = _get_f64(d, "FREFF",               0.07),
        q_re0                   = _get_f64(d, "Q_RE0",               1.0),
        q_v1                    = _get_f64(d, "Q_V1",                0.0),
        q_v2                    = _get_f64(d, "Q_V2",                0.0),
        q_fz2                   = _get_f64(d, "Q_FZ2",               1.0e-4),
        q_fcx                   = _get_f64(d, "Q_FCX",               0.0),
        q_fcy                   = _get_f64(d, "Q_FCY",               0.0),
        q_cam                   = _get_f64(d, "Q_CAM",               0.0),
        pfz1                    = _get_f64(d, "PFZ1",                0.0),
        bottom_offst            = _get_f64(d, "BOTTOM_OFFST",        0.01),
        bottom_stiff            = _get_f64(d, "BOTTOM_STIFF",        2_000_000.0),
        # [STRUCTURAL]
        longitudinal_stiffness  = _get_f64(d, "LONGITUDINAL_STIFFNESS", 300_000.0),
        lateral_stiffness       = _get_f64(d, "LATERAL_STIFFNESS",      100_000.0),
        yaw_stiffness           = _get_f64(d, "YAW_STIFFNESS",          5_000.0),
        freq_long               = _get_f64(d, "FREQ_LONG",              80.0),
        freq_lat                = _get_f64(d, "FREQ_LAT",               40.0),
        freq_yaw                = _get_f64(d, "FREQ_YAW",               50.0),
        freq_windup             = _get_f64(d, "FREQ_WINDUP",            70.0),
        damp_long               = _get_f64(d, "DAMP_LONG",              0.04),
        damp_lat                = _get_f64(d, "DAMP_LAT",               0.04),
        damp_yaw                = _get_f64(d, "DAMP_YAW",               0.04),
        damp_windup             = _get_f64(d, "DAMP_WINDUP",            0.04),
        damp_residual           = _get_f64(d, "DAMP_RESIDUAL",          0.002),
        damp_vlow               = _get_f64(d, "DAMP_VLOW",              0.001),
        q_bvx                   = _get_f64(d, "Q_BVX",                  0.0),
        q_bvt                   = _get_f64(d, "Q_BVT",                  0.0),
        # [TYRE_CARCASS]
        pcfx1                   = _get_f64(d, "PCFX1",               0.0),
        pcfx2                   = _get_f64(d, "PCFX2",               0.0),
        pcfx3                   = _get_f64(d, "PCFX3",               0.0),
        pcfy1                   = _get_f64(d, "PCFY1",               0.0),
        pcfy2                   = _get_f64(d, "PCFY2",               0.0),
        pcfy3                   = _get_f64(d, "PCFY3",               0.0),
        pcmz1                   = _get_f64(d, "PCMZ1",               0.0),
        # [CONTACT_PATCH]
        q_ra1                   = _get_f64(d, "Q_RA1",               0.5),
        q_ra2                   = _get_f64(d, "Q_RA2",               1.0),
        q_rb1                   = _get_f64(d, "Q_RB1",               1.0),
        q_rb2                   = _get_f64(d, "Q_RB2",               -1.0),
        ellips_shift            = _get_f64(d, "ELLIPS_SHIFT",        0.8),
        ellips_length           = _get_f64(d, "ELLIPS_LENGTH",       1.0),
        ellips_height           = _get_f64(d, "ELLIPS_HEIGHT",       1.0),
        ellips_order            = _get_f64(d, "ELLIPS_ORDER",        1.8),
        ellips_max_step         = _get_f64(d, "ELLIPS_MAX_STEP",     0.025),
        ellips_nwidth           = _get_f64(d, "ELLIPS_NWIDTH",       10.0),
        ellips_nlength          = _get_f64(d, "ELLIPS_NLENGTH",      10.0),
        # Ranges
        presmin                 = _get_f64(d, "PRESMIN",             10_000.0),
        presmax                 = _get_f64(d, "PRESMAX",             1_000_000.0),
        fzmin                   = _get_f64(d, "FZMIN",               100.0),
        fzmax                   = _get_f64(d, "FZMAX",               10_000.0),
        kpumin                  = _get_f64(d, "KPUMIN",              -1.5),
        kpumax                  = _get_f64(d, "KPUMAX",               1.5),
        alpmin                  = _get_f64(d, "ALPMIN",              -1.5),
        alpmax                  = _get_f64(d, "ALPMAX",               1.5),
        cammin                  = _get_f64(d, "CAMMIN",              -0.175),
        cammax                  = _get_f64(d, "CAMMAX",               0.175),
        # [SCALING_COEFFICIENTS]
        lfzo                    = _get_f64(d, "LFZO",                1.0),
        lcx                     = _get_f64(d, "LCX",                 1.0),
        lmux                    = _get_f64(d, "LMUX",                1.0),
        lex                     = _get_f64(d, "LEX",                 1.0),
        lkx                     = _get_f64(d, "LKX",                 1.0),
        lhx                     = _get_f64(d, "LHX",                 1.0),
        lvx                     = _get_f64(d, "LVX",                 1.0),
        lcy                     = _get_f64(d, "LCY",                 1.0),
        lmuy                    = _get_f64(d, "LMUY",                1.0),
        ley                     = _get_f64(d, "LEY",                 1.0),
        lky                     = _get_f64(d, "LKY",                 1.0),
        lhy                     = _get_f64(d, "LHY",                 1.0),
        lvy                     = _get_f64(d, "LVY",                 1.0),
        ltr                     = _get_f64(d, "LTR",                 1.0),
        lres                    = _get_f64(d, "LRES",                1.0),
        lxal                    = _get_f64(d, "LXAL",                1.0),
        lyka                    = _get_f64(d, "LYKA",                1.0),
        lvyka                   = _get_f64(d, "LVYKA",               1.0),
        ls                      = _get_f64(d, "LS",                  1.0),
        lkyc                    = _get_f64(d, "LKYC",                1.0),
        lkzc                    = _get_f64(d, "LKZC",                1.0),
        lvmx                    = _get_f64(d, "LVMX",                1.0),
        lmx                     = _get_f64(d, "LMX",                 1.0),
        lmy                     = _get_f64(d, "LMY",                 1.0),
        lmp                     = _get_f64(d, "LMP",                 1.0),
        # [LONGITUDINAL_COEFFICIENTS]
        pcx1                    = _get_f64(d, "PCX1",                1.65),
        pdx1                    = _get_f64(d, "PDX1",                1.3),
        pdx2                    = _get_f64(d, "PDX2",               -0.15),
        pdx3                    = _get_f64(d, "PDX3",                0.0),
        pex1                    = _get_f64(d, "PEX1",                0.0),
        pex2                    = _get_f64(d, "PEX2",                0.0),
        pex3                    = _get_f64(d, "PEX3",                0.0),
        pex4                    = _get_f64(d, "PEX4",                0.0),
        pkx1                    = _get_f64(d, "PKX1",               20.0),
        pkx2                    = _get_f64(d, "PKX2",                0.0),
        pkx3                    = _get_f64(d, "PKX3",                0.0),
        phx1                    = _get_f64(d, "PHX1",                0.0),
        phx2                    = _get_f64(d, "PHX2",                0.0),
        pvx1                    = _get_f64(d, "PVX1",                0.0),
        pvx2                    = _get_f64(d, "PVX2",                0.0),
        ppx1                    = _get_f64(d, "PPX1",                0.0),
        ppx2                    = _get_f64(d, "PPX2",                0.0),
        ppx3                    = _get_f64(d, "PPX3",                0.0),
        ppx4                    = _get_f64(d, "PPX4",                0.0),
        rbx1                    = _get_f64(d, "RBX1",               20.0),
        rbx2                    = _get_f64(d, "RBX2",               15.0),
        rbx3                    = _get_f64(d, "RBX3",                0.0),
        rcx1                    = _get_f64(d, "RCX1",                1.0),
        rex1                    = _get_f64(d, "REX1",                0.0),
        rex2                    = _get_f64(d, "REX2",                0.0),
        rhx1                    = _get_f64(d, "RHX1",                0.0),
        # [OVERTURNING_COEFFICIENTS]
        qsx1                    = _get_f64(d, "QSX1",                0.0),
        qsx2                    = _get_f64(d, "QSX2",                0.0),
        qsx3                    = _get_f64(d, "QSX3",                0.0),
        qsx4                    = _get_f64(d, "QSX4",                5.0),
        qsx5                    = _get_f64(d, "QSX5",                1.0),
        qsx6                    = _get_f64(d, "QSX6",               10.0),
        qsx7                    = _get_f64(d, "QSX7",                0.0),
        qsx8                    = _get_f64(d, "QSX8",                0.0),
        qsx9                    = _get_f64(d, "QSX9",                0.4),
        qsx10                   = _get_f64(d, "QSX10",               0.0),
        qsx11                   = _get_f64(d, "QSX11",               5.0),
        qsx12                   = _get_f64(d, "QSX12",               0.0),
        qsx13                   = _get_f64(d, "QSX13",               0.0),
        qsx14                   = _get_f64(d, "QSX14",               0.0),
        ppmx1                   = _get_f64(d, "PPMX1",               0.0),
        # [LATERAL_COEFFICIENTS]
        pcy1                    = _get_f64(d, "PCY1",                1.3),
        pdy1                    = _get_f64(d, "PDY1",                1.1),
        pdy2                    = _get_f64(d, "PDY2",               -0.15),
        pdy3                    = _get_f64(d, "PDY3",                0.0),
        pey1                    = _get_f64(d, "PEY1",                0.0),
        pey2                    = _get_f64(d, "PEY2",                0.0),
        pey3                    = _get_f64(d, "PEY3",                0.0),
        pey4                    = _get_f64(d, "PEY4",                0.0),
        pey5                    = _get_f64(d, "PEY5",                0.0),
        pky1                    = _get_f64(d, "PKY1",              -20.0),
        pky2                    = _get_f64(d, "PKY2",                1.0),
        pky3                    = _get_f64(d, "PKY3",                0.0),
        pky4                    = _get_f64(d, "PKY4",                2.0),
        pky5                    = _get_f64(d, "PKY5",                0.0),
        pky6                    = _get_f64(d, "PKY6",               -1.0),
        pky7                    = _get_f64(d, "PKY7",                0.0),
        phy1                    = _get_f64(d, "PHY1",                0.0),
        phy2                    = _get_f64(d, "PHY2",                0.0),
        pvy1                    = _get_f64(d, "PVY1",                0.0),
        pvy2                    = _get_f64(d, "PVY2",                0.0),
        pvy3                    = _get_f64(d, "PVY3",                0.0),
        pvy4                    = _get_f64(d, "PVY4",                0.0),
        ppy1                    = _get_f64(d, "PPY1",                0.0),
        ppy2                    = _get_f64(d, "PPY2",                0.0),
        ppy3                    = _get_f64(d, "PPY3",                0.0),
        ppy4                    = _get_f64(d, "PPY4",                0.0),
        ppy5                    = _get_f64(d, "PPY5",                0.0),
        rby1                    = _get_f64(d, "RBY1",               10.0),
        rby2                    = _get_f64(d, "RBY2",               10.0),
        rby3                    = _get_f64(d, "RBY3",                0.0),
        rby4                    = _get_f64(d, "RBY4",                0.0),
        rcy1                    = _get_f64(d, "RCY1",                1.0),
        rey1                    = _get_f64(d, "REY1",                0.0),
        rey2                    = _get_f64(d, "REY2",                0.0),
        rhy1                    = _get_f64(d, "RHY1",                0.0),
        rhy2                    = _get_f64(d, "RHY2",                0.0),
        rvy1                    = _get_f64(d, "RVY1",                0.0),
        rvy2                    = _get_f64(d, "RVY2",                0.0),
        rvy3                    = _get_f64(d, "RVY3",                0.0),
        rvy4                    = _get_f64(d, "RVY4",               20.0),
        rvy5                    = _get_f64(d, "RVY5",                2.0),
        rvy6                    = _get_f64(d, "RVY6",               10.0),
        # [ROLLING_COEFFICIENTS]
        qsy1                    = _get_f64(d, "QSY1",                0.01),
        qsy2                    = _get_f64(d, "QSY2",                0.0),
        qsy3                    = _get_f64(d, "QSY3",                4.0e-4),
        qsy4                    = _get_f64(d, "QSY4",                4.0e-5),
        qsy5                    = _get_f64(d, "QSY5",                0.0),
        qsy6                    = _get_f64(d, "QSY6",                0.0),
        qsy7                    = _get_f64(d, "QSY7",                0.85),
        qsy8                    = _get_f64(d, "QSY8",               -0.4),
        # [ALIGNING_COEFFICIENTS]
        qbz1                    = _get_f64(d, "QBZ1",               10.0),
        qbz2                    = _get_f64(d, "QBZ2",                0.0),
        qbz3                    = _get_f64(d, "QBZ3",                0.0),
        qbz4                    = _get_f64(d, "QBZ4",                0.0),
        qbz5                    = _get_f64(d, "QBZ5",                0.0),
        qbz9                    = _get_f64(d, "QBZ9",               10.0),
        qbz10                   = _get_f64(d, "QBZ10",               0.0),
        qcz1                    = _get_f64(d, "QCZ1",                1.1),
        qdz1                    = _get_f64(d, "QDZ1",                0.12),
        qdz2                    = _get_f64(d, "QDZ2",                0.0),
        qdz3                    = _get_f64(d, "QDZ3",                0.0),
        qdz4                    = _get_f64(d, "QDZ4",                0.0),
        qdz6                    = _get_f64(d, "QDZ6",                0.0),
        qdz7                    = _get_f64(d, "QDZ7",                0.0),
        qdz8                    = _get_f64(d, "QDZ8",               -0.05),
        qdz9                    = _get_f64(d, "QDZ9",                0.0),
        qdz10                   = _get_f64(d, "QDZ10",               0.0),
        qdz11                   = _get_f64(d, "QDZ11",               0.0),
        qez1                    = _get_f64(d, "QEZ1",                0.0),
        qez2                    = _get_f64(d, "QEZ2",                0.0),
        qez3                    = _get_f64(d, "QEZ3",                0.0),
        qez4                    = _get_f64(d, "QEZ4",                0.0),
        qez5                    = _get_f64(d, "QEZ5",                0.0),
        qhz1                    = _get_f64(d, "QHZ1",                0.0),
        qhz2                    = _get_f64(d, "QHZ2",                0.0),
        qhz3                    = _get_f64(d, "QHZ3",                0.0),
        qhz4                    = _get_f64(d, "QHZ4",                0.0),
        ppz1                    = _get_f64(d, "PPZ1",                0.0),
        ppz2                    = _get_f64(d, "PPZ2",                0.0),
        ssz1                    = _get_f64(d, "SSZ1",                0.0),
        ssz2                    = _get_f64(d, "SSZ2",                0.0),
        ssz3                    = _get_f64(d, "SSZ3",                0.0),
        ssz4                    = _get_f64(d, "SSZ4",                0.0),
        # [TURNSLIP_COEFFICIENTS]
        pdxp1                   = _get_f64(d, "PDXP1",               0.4),
        pdxp2                   = _get_f64(d, "PDXP2",               0.0),
        pdxp3                   = _get_f64(d, "PDXP3",               0.0),
        pkyp1                   = _get_f64(d, "PKYP1",               1.0),
        pdyp1                   = _get_f64(d, "PDYP1",               0.4),
        pdyp2                   = _get_f64(d, "PDYP2",               0.0),
        pdyp3                   = _get_f64(d, "PDYP3",               0.0),
        pdyp4                   = _get_f64(d, "PDYP4",               0.0),
        phyp1                   = _get_f64(d, "PHYP1",               1.0),
        phyp2                   = _get_f64(d, "PHYP2",               0.15),
        phyp3                   = _get_f64(d, "PHYP3",               0.0),
        phyp4                   = _get_f64(d, "PHYP4",              -4.0),
        pecp1                   = _get_f64(d, "PECP1",               0.5),
        pecp2                   = _get_f64(d, "PECP2",               0.0),
        qdtp1                   = _get_f64(d, "QDTP1",              10.0),
        qcrp1                   = _get_f64(d, "QCRP1",               0.2),
        qcrp2                   = _get_f64(d, "QCRP2",               0.1),
        qbrp1                   = _get_f64(d, "QBRP1",               0.1),
        qdrp1                   = _get_f64(d, "QDRP1",               1.0),
        # MF 5.2 only
        q_a1                    = _get_f64(d, "Q_A1",                0.0),
        q_a2                    = _get_f64(d, "Q_A2",                0.0),
        phy3                    = _get_f64(d, "PHY3",                0.0),
        ptx1                    = _get_f64(d, "PTX1",                0.0),
        ptx2                    = _get_f64(d, "PTX2",                0.0),
        ptx3                    = _get_f64(d, "PTX3",                0.0),
        pty1                    = _get_f64(d, "PTY1",                0.0),
        pty2                    = _get_f64(d, "PTY2",                0.0),
        lsgkp                   = _get_f64(d, "LSGKP",               1.0),
        lsgal                   = _get_f64(d, "LSGAL",               1.0),
        # MF 6.2 only
        switch_integ            = _get_f64(d, "SWITCH_INTEG",         0.0),
        q_fcy2                  = _get_f64(d, "Q_FCY2",               0.0),
        q_cam1                  = _get_f64(d, "Q_CAM1",               0.0),
        q_cam2                  = _get_f64(d, "Q_CAM2",               0.0),
        q_cam3                  = _get_f64(d, "Q_CAM3",               0.0),
        q_fys1                  = _get_f64(d, "Q_FYS1",               0.0),
        q_fys2                  = _get_f64(d, "Q_FYS2",               0.0),
        q_fys3                  = _get_f64(d, "Q_FYS3",               0.0),
        env_c1                  = _get_f64(d, "ENV_C1",               0.0),
        env_c2                  = _get_f64(d, "ENV_C2",               0.0),
        # Simulation config
        hmax_local              = _get_f64(d, "HMAX_LOCAL",           2.5e-4),
        time_switch_integ       = _get_f64(d, "TIME_SWITCH_INTEG",    0.1),
    )
end

# ------------------------------------------------------------------------------
# Low-level helpers
# ------------------------------------------------------------------------------

@inline function _get_f64(d::Dict{String,Any}, key::String, default::Float64) ::Float64
    haskey(d, key) ? Float64(d[key]) : default
end

@inline function _get_str(d::Dict{String,Any}, key::String, default::String) ::String
    haskey(d, key) ? String(d[key]) : default
end

# Try to parse a string as Float64; if it fails, return it as a String
function _try_parse_float_or_string(s::AbstractString, key::AbstractString) ::Any
    v = tryparse(Float64, s)
    if isnothing(v)
        return _strip_quotes(s)
    end
    return v
end

@inline function _try_parse_float(s::AbstractString, key::AbstractString) ::Float64
    v = tryparse(Float64, s)
    isnothing(v) && error("read_tir: could not parse value for key \"$key\": \"$s\"")
    return v
end

@inline _strip_quotes(s::AbstractString) = strip(s, ['\'', '"', ' '])