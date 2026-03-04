# ==============================================================================
# generate_reference_data.jl
# Run from MFeval.jl/ root:  julia --project=. test/generate_reference_data.jl
# ==============================================================================

include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
using .MFeval

let
    OUTDIR  = joinpath(@__DIR__, "reference_data")
    FIX_DIR = joinpath(@__DIR__, "fixtures")
    COLS    = ["Fx","Fy","Fz","Mx","My","Mz",
               "kappa","alpha","gamma","phit","Vx","pressure",
               "Re","rho","two_a","t","mux","muy","omega","Rl","two_b","Mzr",
               "Cx","Cy","Cz","Kya","sigmax","sigmay","inst_Kya","Kxk"]

    mkpath(OUTDIR)

    function write_csv(path, mat)
        open(path, "w") do io
            println(io, join(COLS, ","))
            for i in 1:size(mat, 1)
                println(io, join(mat[i, :], ","))
            end
        end
        println("  wrote $(basename(path))  ($(size(mat,1)) rows)")
    end

    p61 = read_tir(joinpath(FIX_DIR, "MagicFormula61_Parameters.tir"))
    p62 = read_tir(joinpath(FIX_DIR, "MagicFormula62_Parameters.tir"))
    p52 = read_tir(joinpath(FIX_DIR, "MagicFormula52_Parameters.tir"))
    M   = MFModes(111)

    Fnom61 = p61.fnomin;  Vnom61 = p61.longvl
    Fnom62 = p62.fnomin;  Vnom62 = p62.longvl
    Fnom52 = p52.fnomin;  Vnom52 = p52.longvl

    println("Generating reference data…")

    # helpers
    rep(x, n) = fill(Float64(x), n)
    inp6(Fz, κ, α, γ, Vx) = [Fz κ α γ zeros(length(Fz)) Vx]

    # MF6.1 — pure lateral
    let N=200, α=collect(range(-0.5,0.5,length=200))
        write_csv(joinpath(OUTDIR,"mf61_pure_lat.csv"),
            mfeval(p61, inp6(rep(Fnom61,N),zeros(N),α,zeros(N),rep(Vnom61,N)), M))
    end

    # MF6.1 — pure longitudinal
    let N=200, κ=collect(range(-0.5,0.5,length=200))
        write_csv(joinpath(OUTDIR,"mf61_pure_lon.csv"),
            mfeval(p61, inp6(rep(Fnom61,N),κ,zeros(N),zeros(N),rep(Vnom61,N)), M))
    end

    # MF6.1 — combined slip (11×11)
    let κs=collect(range(-0.3,0.3,length=11)), αs=collect(range(-0.3,0.3,length=11))
        κv=[κ for κ in κs for _ in αs]; αv=[α for _ in κs for α in αs]; N=length(κv)
        write_csv(joinpath(OUTDIR,"mf61_combined.csv"),
            mfeval(p61, inp6(rep(Fnom61,N),κv,αv,zeros(N),rep(Vnom61,N)), M))
    end

    # MF6.1 — camber sweep
    let N=50, γ=collect(range(-p61.cammax,p61.cammax,length=50))
        write_csv(joinpath(OUTDIR,"mf61_camber.csv"),
            mfeval(p61, inp6(rep(Fnom61,N),zeros(N),zeros(N),γ,rep(Vnom61,N)), M))
    end

    # MF6.1 — load sweep
    let N=50, Fz=collect(range(0.25*Fnom61,2.0*Fnom61,length=50))
        write_csv(joinpath(OUTDIR,"mf61_load.csv"),
            mfeval(p61, inp6(Fz,zeros(N),rep(0.1,N),zeros(N),rep(Vnom61,N)), M))
    end

    # MF6.1 — pressure sweep (7-column input)
    let N=30, pres=collect(range(p61.presmin*1.1,p61.presmax*0.9,length=30))
        mat7 = hcat(rep(Fnom61,N),zeros(N),rep(0.1,N),zeros(N),zeros(N),rep(Vnom61,N),pres)
        write_csv(joinpath(OUTDIR,"mf61_pressure.csv"), mfeval(p61, mat7, M))
    end

    # MF6.2 — pure lateral
    let N=200, α=collect(range(-0.5,0.5,length=200))
        write_csv(joinpath(OUTDIR,"mf62_pure_lat.csv"),
            mfeval(p62, inp6(rep(Fnom62,N),zeros(N),α,zeros(N),rep(Vnom62,N)), M))
    end

    # MF5.2 — pure lateral
    let N=200, α=collect(range(-0.5,0.5,length=200))
        write_csv(joinpath(OUTDIR,"mf52_pure_lat.csv"),
            mfeval(p52, inp6(rep(Fnom52,N),zeros(N),α,zeros(N),rep(Vnom52,N)), M))
    end

    # MF6.1 — edge cases
    let e = Float64[
            0.0         0.0              0.0          0.0          0.0  Vnom61 ;
            Fnom61      0.0              0.0          0.0          0.0  0.0    ;
            Fnom61     -1.0              0.0          0.0          0.0  Vnom61 ;
            Fnom61      0.0              0.0          p61.cammax   0.0  Vnom61 ;
            Fnom61      0.0              0.0         -p61.cammax   0.0  Vnom61 ;
            Fnom61      0.0              p61.alpmax*1.5  0.0       0.0  Vnom61 ;
            Fnom61      p61.kpumax*1.5   0.0          0.0          0.0  Vnom61 ;
            Fnom61      0.0              0.1          0.0          0.0  0.3    ;
            Fnom61*2.0  0.0              0.1          0.0          0.0  Vnom61 ;
        ]
        write_csv(joinpath(OUTDIR,"mf61_edge_cases.csv"), mfeval(p61, e, M))
    end

    println("Done. Reference CSVs written to test/reference_data/")
end