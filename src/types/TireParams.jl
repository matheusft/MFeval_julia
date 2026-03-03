# ==============================================================================
# TireParams.jl
#
# Immutable, fully-concrete struct holding all Magic Formula tyre coefficients.
#
# The type parameter V encodes the MF version at compile time:
#   V = :MF52  → FITTYP 6 or 21
#   V = :MF61  → FITTYP 61
#   V = :MF62  → FITTYP 62
#
# Every field is Float64. String metadata (TYRESIDE, etc.) is stored separately
# in TireMetadata to keep the hot-path struct free of non-numeric data.
#
# Design note: all fields have default values matching the readTIR.m defaults
# so that a TireParams can be constructed incrementally from a Dict during
# TIR parsing without requiring all keys to be present in every file.
# ==============================================================================

# ------------------------------------------------------------------------------
# Metadata — strings and integers that are never used in arithmetic
# ------------------------------------------------------------------------------
struct TireMetadata
    tyreside             ::String   # 'LEFT' or 'RIGHT'
    file_type            ::String
    file_version         ::Float64
    file_format          ::String
    property_file_format ::String
    function_name        ::String   # MF6.2 only
    units_length         ::String   # e.g. 'meter'
    units_force          ::String   # e.g. 'newton'
    units_angle          ::String   # e.g. 'radians'
    fittyp               ::Int      # Raw FITTYP integer from TIR file
    use_mode             ::Int
    n_tire_states        ::Int
    user_sub_id          ::Int
end

function TireMetadata(;
    tyreside             = "LEFT",
    file_type            = "default",
    file_version         = 3.0,
    file_format          = "default",
    property_file_format = "default",
    function_name        = "default",
    units_length         = "meter",
    units_force          = "newton",
    units_angle          = "radians",
    fittyp               = 61,
    use_mode             = 124,
    n_tire_states        = 4,
    user_sub_id          = 815,
)
    TireMetadata(tyreside, file_type, file_version, file_format,
                 property_file_format, function_name,
                 units_length, units_force, units_angle,
                 fittyp, use_mode, n_tire_states, user_sub_id)
end

# ------------------------------------------------------------------------------
# TireParams{V} — all numeric coefficients, grouped by TIR section
# ------------------------------------------------------------------------------
struct TireParams{V}

    metadata ::TireMetadata

    # ── [MODEL] ────────────────────────────────────────────────────────────────
    longvl           ::Float64   # Nominal speed                         [m/s]
    vxlow            ::Float64   # Lower boundary of slip calculation    [m/s]
    road_increment   ::Float64
    road_direction   ::Float64

    # ── [DIMENSION] ───────────────────────────────────────────────────────────
    unloaded_radius  ::Float64   # R0: free tyre radius                  [m]
    width            ::Float64   # Nominal tyre width                    [m]
    aspect_ratio     ::Float64
    rim_radius       ::Float64   # Rrim                                  [m]
    rim_width        ::Float64   # [m]

    # ── [OPERATING_CONDITIONS] ────────────────────────────────────────────────
    inflpres         ::Float64   # Inflation pressure (operating)        [Pa]
    nompres          ::Float64   # Nominal inflation pressure pi0        [Pa]

    # ── [INERTIA] ─────────────────────────────────────────────────────────────
    mass             ::Float64   # [kg]
    ixx              ::Float64   # [kg·m²]
    iyy              ::Float64   # [kg·m²]
    belt_mass        ::Float64
    belt_ixx         ::Float64
    belt_iyy         ::Float64
    gravity          ::Float64

    # ── [VERTICAL] ────────────────────────────────────────────────────────────
    fnomin           ::Float64   # Fz0: nominal wheel load               [N]
    vertical_stiffness ::Float64 # Cz0                                   [N/m]
    vertical_damping   ::Float64
    mc_contour_a     ::Float64
    mc_contour_b     ::Float64
    breff            ::Float64   # Effective rolling radius parameter
    dreff            ::Float64
    freff            ::Float64
    q_re0            ::Float64
    q_v1             ::Float64
    q_v2             ::Float64
    q_fz2            ::Float64
    q_fcx            ::Float64
    q_fcy            ::Float64
    q_cam            ::Float64
    pfz1             ::Float64   # Pressure effect on vertical stiffness
    bottom_offst     ::Float64
    bottom_stiff     ::Float64

    # ── [STRUCTURAL] ──────────────────────────────────────────────────────────
    longitudinal_stiffness ::Float64  # cx0
    lateral_stiffness      ::Float64  # cy0
    yaw_stiffness          ::Float64
    freq_long              ::Float64
    freq_lat               ::Float64
    freq_yaw               ::Float64
    freq_windup            ::Float64
    damp_long              ::Float64
    damp_lat               ::Float64
    damp_yaw               ::Float64
    damp_windup            ::Float64
    damp_residual          ::Float64
    damp_vlow              ::Float64
    q_bvx                  ::Float64
    q_bvt                  ::Float64

    # ── [TYRE_CARCASS] (relaxation stiffnesses) ───────────────────────────────
    pcfx1            ::Float64
    pcfx2            ::Float64
    pcfx3            ::Float64
    pcfy1            ::Float64
    pcfy2            ::Float64
    pcfy3            ::Float64
    pcmz1            ::Float64

    # ── [CONTACT_PATCH] ───────────────────────────────────────────────────────
    q_ra1            ::Float64
    q_ra2            ::Float64
    q_rb1            ::Float64
    q_rb2            ::Float64
    ellips_shift     ::Float64
    ellips_length    ::Float64
    ellips_height    ::Float64
    ellips_order     ::Float64
    ellips_max_step  ::Float64
    ellips_nwidth    ::Float64
    ellips_nlength   ::Float64

    # ── [INFLATION_PRESSURE_RANGE] ────────────────────────────────────────────
    presmin          ::Float64
    presmax          ::Float64

    # ── [VERTICAL_FORCE_RANGE] ────────────────────────────────────────────────
    fzmin            ::Float64
    fzmax            ::Float64

    # ── [LONG_SLIP_RANGE] ─────────────────────────────────────────────────────
    kpumin           ::Float64
    kpumax           ::Float64

    # ── [SLIP_ANGLE_RANGE] ────────────────────────────────────────────────────
    alpmin           ::Float64
    alpmax           ::Float64

    # ── [INCLINATION_ANGLE_RANGE] ─────────────────────────────────────────────
    cammin           ::Float64
    cammax           ::Float64

    # ── [SCALING_COEFFICIENTS] ────────────────────────────────────────────────
    lfzo             ::Float64
    lcx              ::Float64
    lmux             ::Float64
    lex              ::Float64
    lkx              ::Float64
    lhx              ::Float64
    lvx              ::Float64
    lcy              ::Float64
    lmuy             ::Float64
    ley              ::Float64
    lky              ::Float64
    lhy              ::Float64
    lvy              ::Float64
    ltr              ::Float64
    lres             ::Float64
    lxal             ::Float64
    lyka             ::Float64
    lvyka            ::Float64
    ls               ::Float64
    lkyc             ::Float64
    lkzc             ::Float64
    lvmx             ::Float64
    lmx              ::Float64
    lmy              ::Float64
    lmp              ::Float64

    # ── [LONGITUDINAL_COEFFICIENTS] ───────────────────────────────────────────
    pcx1             ::Float64
    pdx1             ::Float64
    pdx2             ::Float64
    pdx3             ::Float64
    pex1             ::Float64
    pex2             ::Float64
    pex3             ::Float64
    pex4             ::Float64
    pkx1             ::Float64
    pkx2             ::Float64
    pkx3             ::Float64
    phx1             ::Float64
    phx2             ::Float64
    pvx1             ::Float64
    pvx2             ::Float64
    ppx1             ::Float64
    ppx2             ::Float64
    ppx3             ::Float64
    ppx4             ::Float64
    rbx1             ::Float64
    rbx2             ::Float64
    rbx3             ::Float64
    rcx1             ::Float64
    rex1             ::Float64
    rex2             ::Float64
    rhx1             ::Float64

    # ── [OVERTURNING_COEFFICIENTS] ────────────────────────────────────────────
    qsx1             ::Float64
    qsx2             ::Float64
    qsx3             ::Float64
    qsx4             ::Float64
    qsx5             ::Float64
    qsx6             ::Float64
    qsx7             ::Float64
    qsx8             ::Float64
    qsx9             ::Float64
    qsx10            ::Float64
    qsx11            ::Float64
    qsx12            ::Float64
    qsx13            ::Float64
    qsx14            ::Float64
    ppmx1            ::Float64

    # ── [LATERAL_COEFFICIENTS] ────────────────────────────────────────────────
    pcy1             ::Float64
    pdy1             ::Float64
    pdy2             ::Float64
    pdy3             ::Float64
    pey1             ::Float64
    pey2             ::Float64
    pey3             ::Float64
    pey4             ::Float64
    pey5             ::Float64
    pky1             ::Float64
    pky2             ::Float64
    pky3             ::Float64
    pky4             ::Float64
    pky5             ::Float64
    pky6             ::Float64
    pky7             ::Float64
    phy1             ::Float64
    phy2             ::Float64
    pvy1             ::Float64
    pvy2             ::Float64
    pvy3             ::Float64
    pvy4             ::Float64
    ppy1             ::Float64
    ppy2             ::Float64
    ppy3             ::Float64
    ppy4             ::Float64
    ppy5             ::Float64
    rby1             ::Float64
    rby2             ::Float64
    rby3             ::Float64
    rby4             ::Float64
    rcy1             ::Float64
    rey1             ::Float64
    rey2             ::Float64
    rhy1             ::Float64
    rhy2             ::Float64
    rvy1             ::Float64
    rvy2             ::Float64
    rvy3             ::Float64
    rvy4             ::Float64
    rvy5             ::Float64
    rvy6             ::Float64

    # ── [ROLLING_COEFFICIENTS] ────────────────────────────────────────────────
    qsy1             ::Float64
    qsy2             ::Float64
    qsy3             ::Float64
    qsy4             ::Float64
    qsy5             ::Float64
    qsy6             ::Float64
    qsy7             ::Float64
    qsy8             ::Float64

    # ── [ALIGNING_COEFFICIENTS] ───────────────────────────────────────────────
    qbz1             ::Float64
    qbz2             ::Float64
    qbz3             ::Float64
    qbz4             ::Float64
    qbz5             ::Float64
    qbz9             ::Float64
    qbz10            ::Float64
    qcz1             ::Float64
    qdz1             ::Float64
    qdz2             ::Float64
    qdz3             ::Float64
    qdz4             ::Float64
    qdz6             ::Float64
    qdz7             ::Float64
    qdz8             ::Float64
    qdz9             ::Float64
    qdz10            ::Float64
    qdz11            ::Float64
    qez1             ::Float64
    qez2             ::Float64
    qez3             ::Float64
    qez4             ::Float64
    qez5             ::Float64
    qhz1             ::Float64
    qhz2             ::Float64
    qhz3             ::Float64
    qhz4             ::Float64
    ppz1             ::Float64
    ppz2             ::Float64
    ssz1             ::Float64
    ssz2             ::Float64
    ssz3             ::Float64
    ssz4             ::Float64

    # ── [TURNSLIP_COEFFICIENTS] ───────────────────────────────────────────────
    pdxp1            ::Float64
    pdxp2            ::Float64
    pdxp3            ::Float64
    pkyp1            ::Float64
    pdyp1            ::Float64
    pdyp2            ::Float64
    pdyp3            ::Float64
    pdyp4            ::Float64
    phyp1            ::Float64
    phyp2            ::Float64
    phyp3            ::Float64
    phyp4            ::Float64
    pecp1            ::Float64
    pecp2            ::Float64
    qdtp1            ::Float64
    qcrp1            ::Float64
    qcrp2            ::Float64
    qbrp1            ::Float64
    qdrp1            ::Float64

    # ── MF 5.2 only ───────────────────────────────────────────────────────────
    q_a1             ::Float64
    q_a2             ::Float64
    phy3             ::Float64
    ptx1             ::Float64
    ptx2             ::Float64
    ptx3             ::Float64
    pty1             ::Float64
    pty2             ::Float64
    lsgkp            ::Float64
    lsgal            ::Float64

    # ── MF 6.2 only ───────────────────────────────────────────────────────────
    switch_integ     ::Float64
    q_fcy2           ::Float64
    q_cam1           ::Float64
    q_cam2           ::Float64
    q_cam3           ::Float64
    q_fys1           ::Float64
    q_fys2           ::Float64
    q_fys3           ::Float64
    env_c1           ::Float64
    env_c2           ::Float64

    # ── Simulation config ─────────────────────────────────────────────────────
    hmax_local       ::Float64
    time_switch_integ ::Float64
end

# ------------------------------------------------------------------------------
# Helper: resolve version symbol from raw FITTYP integer
# ------------------------------------------------------------------------------
function fittyp_to_version(fittyp::Int) ::Symbol
    if fittyp == 6 || fittyp == 21
        return :MF52
    elseif fittyp == 61
        return :MF61
    elseif fittyp == 62
        return :MF62
    else
        error("TireParams: unsupported FITTYP = $fittyp. " *
              "Expected 6, 21 (MF5.2), 61 (MF6.1), or 62 (MF6.2).")
    end
end

# ------------------------------------------------------------------------------
# Keyword constructor with all defaults matching readTIR.m
# The caller (read_tir) passes only the keys present in the file; everything
# else falls back to the MATLAB-equivalent defaults.
# ------------------------------------------------------------------------------
function TireParams(;
    metadata                  = TireMetadata(),
    # [MODEL]
    longvl                    = 16.7,
    vxlow                     = 1.0,
    road_increment            = 0.01,
    road_direction            = 1.0,
    # [DIMENSION]
    unloaded_radius           = 0.3135,
    width                     = 0.205,
    aspect_ratio              = 0.6,
    rim_radius                = 0.1905,
    rim_width                 = 0.152,
    # [OPERATING_CONDITIONS]
    inflpres                  = 220000.0,
    nompres                   = 220000.0,
    # [INERTIA]
    mass                      = 9.3,
    ixx                       = 0.391,
    iyy                       = 0.736,
    belt_mass                 = 7.0,
    belt_ixx                  = 0.34,
    belt_iyy                  = 0.6,
    gravity                   = -9.81,
    # [VERTICAL]
    fnomin                    = 4000.0,
    vertical_stiffness        = 200000.0,
    vertical_damping          = 50.0,
    mc_contour_a              = 0.5,
    mc_contour_b              = 0.5,
    breff                     = 8.4,
    dreff                     = 0.27,
    freff                     = 0.07,
    q_re0                     = 1.0,
    q_v1                      = 0.0,
    q_v2                      = 0.0,
    q_fz2                     = 1.0e-4,
    q_fcx                     = 0.0,
    q_fcy                     = 0.0,
    q_cam                     = 0.0,
    pfz1                      = 0.0,
    bottom_offst              = 0.01,
    bottom_stiff              = 2_000_000.0,
    # [STRUCTURAL]
    longitudinal_stiffness    = 300_000.0,
    lateral_stiffness         = 100_000.0,
    yaw_stiffness             = 5_000.0,
    freq_long                 = 80.0,
    freq_lat                  = 40.0,
    freq_yaw                  = 50.0,
    freq_windup               = 70.0,
    damp_long                 = 0.04,
    damp_lat                  = 0.04,
    damp_yaw                  = 0.04,
    damp_windup               = 0.04,
    damp_residual             = 0.002,
    damp_vlow                 = 0.001,
    q_bvx                     = 0.0,
    q_bvt                     = 0.0,
    # [TYRE_CARCASS]
    pcfx1                     = 0.0,
    pcfx2                     = 0.0,
    pcfx3                     = 0.0,
    pcfy1                     = 0.0,
    pcfy2                     = 0.0,
    pcfy3                     = 0.0,
    pcmz1                     = 0.0,
    # [CONTACT_PATCH]
    q_ra1                     = 0.5,
    q_ra2                     = 1.0,
    q_rb1                     = 1.0,
    q_rb2                     = -1.0,
    ellips_shift              = 0.8,
    ellips_length             = 1.0,
    ellips_height             = 1.0,
    ellips_order              = 1.8,
    ellips_max_step           = 0.025,
    ellips_nwidth             = 10.0,
    ellips_nlength            = 10.0,
    # [INFLATION_PRESSURE_RANGE]
    presmin                   = 10_000.0,
    presmax                   = 1_000_000.0,
    # [VERTICAL_FORCE_RANGE]
    fzmin                     = 100.0,
    fzmax                     = 10_000.0,
    # [LONG_SLIP_RANGE]
    kpumin                    = -1.5,
    kpumax                    = 1.5,
    # [SLIP_ANGLE_RANGE]
    alpmin                    = -1.5,
    alpmax                    = 1.5,
    # [INCLINATION_ANGLE_RANGE]
    cammin                    = -0.175,
    cammax                    = 0.175,
    # [SCALING_COEFFICIENTS]
    lfzo                      = 1.0,
    lcx                       = 1.0,
    lmux                      = 1.0,
    lex                       = 1.0,
    lkx                       = 1.0,
    lhx                       = 1.0,
    lvx                       = 1.0,
    lcy                       = 1.0,
    lmuy                      = 1.0,
    ley                       = 1.0,
    lky                       = 1.0,
    lhy                       = 1.0,
    lvy                       = 1.0,
    ltr                       = 1.0,
    lres                      = 1.0,
    lxal                      = 1.0,
    lyka                      = 1.0,
    lvyka                     = 1.0,
    ls                        = 1.0,
    lkyc                      = 1.0,
    lkzc                      = 1.0,
    lvmx                      = 1.0,
    lmx                       = 1.0,
    lmy                       = 1.0,
    lmp                       = 1.0,
    # [LONGITUDINAL_COEFFICIENTS]
    pcx1                      = 1.65,
    pdx1                      = 1.3,
    pdx2                      = -0.15,
    pdx3                      = 0.0,
    pex1                      = 0.0,
    pex2                      = 0.0,
    pex3                      = 0.0,
    pex4                      = 0.0,
    pkx1                      = 20.0,
    pkx2                      = 0.0,
    pkx3                      = 0.0,
    phx1                      = 0.0,
    phx2                      = 0.0,
    pvx1                      = 0.0,
    pvx2                      = 0.0,
    ppx1                      = 0.0,
    ppx2                      = 0.0,
    ppx3                      = 0.0,
    ppx4                      = 0.0,
    rbx1                      = 20.0,
    rbx2                      = 15.0,
    rbx3                      = 0.0,
    rcx1                      = 1.0,
    rex1                      = 0.0,
    rex2                      = 0.0,
    rhx1                      = 0.0,
    # [OVERTURNING_COEFFICIENTS]
    qsx1                      = 0.0,
    qsx2                      = 0.0,
    qsx3                      = 0.0,
    qsx4                      = 5.0,
    qsx5                      = 1.0,
    qsx6                      = 10.0,
    qsx7                      = 0.0,
    qsx8                      = 0.0,
    qsx9                      = 0.4,
    qsx10                     = 0.0,
    qsx11                     = 5.0,
    qsx12                     = 0.0,
    qsx13                     = 0.0,
    qsx14                     = 0.0,
    ppmx1                     = 0.0,
    # [LATERAL_COEFFICIENTS]
    pcy1                      = 1.3,
    pdy1                      = 1.1,
    pdy2                      = -0.15,
    pdy3                      = 0.0,
    pey1                      = 0.0,
    pey2                      = 0.0,
    pey3                      = 0.0,
    pey4                      = 0.0,
    pey5                      = 0.0,
    pky1                      = -20.0,
    pky2                      = 1.0,
    pky3                      = 0.0,
    pky4                      = 2.0,
    pky5                      = 0.0,
    pky6                      = -1.0,
    pky7                      = 0.0,
    phy1                      = 0.0,
    phy2                      = 0.0,
    pvy1                      = 0.0,
    pvy2                      = 0.0,
    pvy3                      = 0.0,
    pvy4                      = 0.0,
    ppy1                      = 0.0,
    ppy2                      = 0.0,
    ppy3                      = 0.0,
    ppy4                      = 0.0,
    ppy5                      = 0.0,
    rby1                      = 10.0,
    rby2                      = 10.0,
    rby3                      = 0.0,
    rby4                      = 0.0,
    rcy1                      = 1.0,
    rey1                      = 0.0,
    rey2                      = 0.0,
    rhy1                      = 0.0,
    rhy2                      = 0.0,
    rvy1                      = 0.0,
    rvy2                      = 0.0,
    rvy3                      = 0.0,
    rvy4                      = 20.0,
    rvy5                      = 2.0,
    rvy6                      = 10.0,
    # [ROLLING_COEFFICIENTS]
    qsy1                      = 0.01,
    qsy2                      = 0.0,
    qsy3                      = 4.0e-4,
    qsy4                      = 4.0e-5,
    qsy5                      = 0.0,
    qsy6                      = 0.0,
    qsy7                      = 0.85,
    qsy8                      = -0.4,
    # [ALIGNING_COEFFICIENTS]
    qbz1                      = 10.0,
    qbz2                      = 0.0,
    qbz3                      = 0.0,
    qbz4                      = 0.0,
    qbz5                      = 0.0,
    qbz9                      = 10.0,
    qbz10                     = 0.0,
    qcz1                      = 1.1,
    qdz1                      = 0.12,
    qdz2                      = 0.0,
    qdz3                      = 0.0,
    qdz4                      = 0.0,
    qdz6                      = 0.0,
    qdz7                      = 0.0,
    qdz8                      = -0.05,
    qdz9                      = 0.0,
    qdz10                     = 0.0,
    qdz11                     = 0.0,
    qez1                      = 0.0,
    qez2                      = 0.0,
    qez3                      = 0.0,
    qez4                      = 0.0,
    qez5                      = 0.0,
    qhz1                      = 0.0,
    qhz2                      = 0.0,
    qhz3                      = 0.0,
    qhz4                      = 0.0,
    ppz1                      = 0.0,
    ppz2                      = 0.0,
    ssz1                      = 0.0,
    ssz2                      = 0.0,
    ssz3                      = 0.0,
    ssz4                      = 0.0,
    # [TURNSLIP_COEFFICIENTS]
    pdxp1                     = 0.4,
    pdxp2                     = 0.0,
    pdxp3                     = 0.0,
    pkyp1                     = 1.0,
    pdyp1                     = 0.4,
    pdyp2                     = 0.0,
    pdyp3                     = 0.0,
    pdyp4                     = 0.0,
    phyp1                     = 1.0,
    phyp2                     = 0.15,
    phyp3                     = 0.0,
    phyp4                     = -4.0,
    pecp1                     = 0.5,
    pecp2                     = 0.0,
    qdtp1                     = 10.0,
    qcrp1                     = 0.2,
    qcrp2                     = 0.1,
    qbrp1                     = 0.1,
    qdrp1                     = 1.0,
    # MF 5.2 only
    q_a1                      = 0.0,
    q_a2                      = 0.0,
    phy3                      = 0.0,
    ptx1                      = 0.0,
    ptx2                      = 0.0,
    ptx3                      = 0.0,
    pty1                      = 0.0,
    pty2                      = 0.0,
    lsgkp                     = 1.0,
    lsgal                     = 1.0,
    # MF 6.2 only
    switch_integ              = 0.0,
    q_fcy2                    = 0.0,
    q_cam1                    = 0.0,
    q_cam2                    = 0.0,
    q_cam3                    = 0.0,
    q_fys1                    = 0.0,
    q_fys2                    = 0.0,
    q_fys3                    = 0.0,
    env_c1                    = 0.0,
    env_c2                    = 0.0,
    # Simulation config
    hmax_local                = 2.5e-4,
    time_switch_integ         = 0.1,
)
    V = fittyp_to_version(metadata.fittyp)
    TireParams{V}(
        metadata,
        longvl, vxlow, road_increment, road_direction,
        unloaded_radius, width, aspect_ratio, rim_radius, rim_width,
        inflpres, nompres,
        mass, ixx, iyy, belt_mass, belt_ixx, belt_iyy, gravity,
        fnomin, vertical_stiffness, vertical_damping,
        mc_contour_a, mc_contour_b,
        breff, dreff, freff,
        q_re0, q_v1, q_v2, q_fz2, q_fcx, q_fcy, q_cam,
        pfz1, bottom_offst, bottom_stiff,
        longitudinal_stiffness, lateral_stiffness, yaw_stiffness,
        freq_long, freq_lat, freq_yaw, freq_windup,
        damp_long, damp_lat, damp_yaw, damp_windup,
        damp_residual, damp_vlow, q_bvx, q_bvt,
        pcfx1, pcfx2, pcfx3, pcfy1, pcfy2, pcfy3, pcmz1,
        q_ra1, q_ra2, q_rb1, q_rb2,
        ellips_shift, ellips_length, ellips_height,
        ellips_order, ellips_max_step, ellips_nwidth, ellips_nlength,
        presmin, presmax,
        fzmin, fzmax,
        kpumin, kpumax,
        alpmin, alpmax,
        cammin, cammax,
        lfzo, lcx, lmux, lex, lkx, lhx, lvx,
        lcy, lmuy, ley, lky, lhy, lvy,
        ltr, lres, lxal, lyka, lvyka, ls, lkyc, lkzc, lvmx, lmx, lmy, lmp,
        pcx1, pdx1, pdx2, pdx3,
        pex1, pex2, pex3, pex4,
        pkx1, pkx2, pkx3,
        phx1, phx2, pvx1, pvx2,
        ppx1, ppx2, ppx3, ppx4,
        rbx1, rbx2, rbx3, rcx1, rex1, rex2, rhx1,
        qsx1, qsx2, qsx3, qsx4, qsx5, qsx6, qsx7,
        qsx8, qsx9, qsx10, qsx11, qsx12, qsx13, qsx14, ppmx1,
        pcy1, pdy1, pdy2, pdy3,
        pey1, pey2, pey3, pey4, pey5,
        pky1, pky2, pky3, pky4, pky5, pky6, pky7,
        phy1, phy2, pvy1, pvy2, pvy3, pvy4,
        ppy1, ppy2, ppy3, ppy4, ppy5,
        rby1, rby2, rby3, rby4, rcy1, rey1, rey2,
        rhy1, rhy2, rvy1, rvy2, rvy3, rvy4, rvy5, rvy6,
        qsy1, qsy2, qsy3, qsy4, qsy5, qsy6, qsy7, qsy8,
        qbz1, qbz2, qbz3, qbz4, qbz5, qbz9, qbz10,
        qcz1,
        qdz1, qdz2, qdz3, qdz4, qdz6, qdz7, qdz8, qdz9, qdz10, qdz11,
        qez1, qez2, qez3, qez4, qez5,
        qhz1, qhz2, qhz3, qhz4,
        ppz1, ppz2,
        ssz1, ssz2, ssz3, ssz4,
        pdxp1, pdxp2, pdxp3, pkyp1,
        pdyp1, pdyp2, pdyp3, pdyp4,
        phyp1, phyp2, phyp3, phyp4,
        pecp1, pecp2,
        qdtp1, qcrp1, qcrp2, qbrp1, qdrp1,
        q_a1, q_a2, phy3, ptx1, ptx2, ptx3, pty1, pty2, lsgkp, lsgal,
        switch_integ, q_fcy2, q_cam1, q_cam2, q_cam3,
        q_fys1, q_fys2, q_fys3, env_c1, env_c2,
        hmax_local, time_switch_integ,
    )
end

# Convenience accessors
version(::TireParams{V}) where {V} = V
is_mf52(p::TireParams) = version(p) === :MF52
is_mf61(p::TireParams) = version(p) === :MF61
is_mf62(p::TireParams) = version(p) === :MF62