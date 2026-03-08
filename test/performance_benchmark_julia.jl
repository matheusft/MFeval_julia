# performance_benchmark_julia.jl
#
# Comprehensive performance benchmarking suite for Julia MFeval.jl
# Generates timing data across various scenarios for comparison with MATLAB
#
# Prerequisites:
#   1. Julia with MFeval.jl package
#   2. TIR files in test/fixtures/
#
# Output:
#   test/julia_performance_results.csv - Detailed timing data
#
# Usage:
#   julia --project=. test/performance_benchmark_julia.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "MFeval.jl"))
using .MFeval
using BenchmarkTools
using Printf
using Statistics
using Dates
using Random

println("="^65)
println("  Julia MFeval.jl Performance Benchmark Suite")
println("  Platform: Julia $(VERSION)")
println("  Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
println("  Threads: $(Threads.nthreads())")
println("="^65)
println()

# Configuration
const FIX_DIR = joinpath(@__DIR__, "fixtures")
const RESULTS_FILE = joinpath(@__DIR__, "julia_performance_results.csv")

# Check if fixtures exist
if !isdir(FIX_DIR)
    error("Fixtures directory not found: $FIX_DIR")
end

# Load TIR files
println("Loading TIR files...")
tir61_path = joinpath(FIX_DIR, "MagicFormula61_Parameters.tir")
tir62_path = joinpath(FIX_DIR, "MagicFormula62_Parameters.tir")
tir52_path = joinpath(FIX_DIR, "MagicFormula52_Parameters.tir")

for path in [tir61_path, tir62_path, tir52_path]
    if !isfile(path)
        error("TIR file not found: $path")
    end
end

const p61 = read_tir(tir61_path)
const p62 = read_tir(tir62_path)
const p52 = read_tir(tir52_path)

@printf("  MF 6.1: Fnom=%.1f N, Vnom=%.1f m/s\n", p61.fnomin, p61.longvl)
@printf("  MF 6.2: Fnom=%.1f N, Vnom=%.1f m/s\n", p62.fnomin, p62.longvl)
@printf("  MF 5.2: Fnom=%.1f N, Vnom=%.1f m/s\n", p52.fnomin, p52.longvl)

# Benchmark configuration
const useMode = 111
const M111 = MFModes(useMode)

# Results storage
results = Vector{NamedTuple{(:test_name, :n_points, :mean_time_us, :min_time_us, :max_time_us, 
                            :std_time_us, :p95_time_us, :category), 
                           Tuple{String, Int, Float64, Float64, Float64, Float64, Float64, String}}}()

# Benchmark helper function
function benchmark_function(f, name::String, n_points::Int, category::String; 
                           evals=1000, samples=1000, seconds=5)
    # Warmup
    for _ in 1:100
        f()
    end
    
    # Main benchmark
    b = @benchmark ($f()) evals=evals samples=samples seconds=seconds
    
    # Extract statistics (convert to microseconds)
    mean_time = median(b).time / 1e3  # ns to μs
    min_time = minimum(b).time / 1e3
    max_time = maximum(b).time / 1e3
    std_time = std(b.times) / 1e3
    p95_time = quantile(b.times, 0.95) / 1e3
    
    @printf("  %-25s (N=%4d): %.2f μs (min=%.2f, max=%.2f, std=%.2f)\n", 
            name, n_points, mean_time, min_time, max_time, std_time)
    
    # Store result
    push!(results, (test_name=name, n_points=n_points, mean_time_us=mean_time, 
                   min_time_us=min_time, max_time_us=max_time, std_time_us=std_time,
                   p95_time_us=p95_time, category=category))
    
    return mean_time
end

## ========================================================================
## BENCHMARK 1: Single Point Evaluation (All MF Versions)  
## ========================================================================
println("\n--- Benchmark 1: Single Point Evaluation ---")

test_cases = [
    ("MF61_single", p61, MFInputs(p61.fnomin, -0.05, 0.1, 0.02, 0.0, p61.longvl)),
    ("MF62_single", p62, MFInputs(p62.fnomin, -0.05, 0.1, 0.02, 0.0, p62.longvl)),
    ("MF52_single", p52, MFInputs(p52.fnomin, -0.05, 0.1, 0.02, 0.0, p52.longvl))
]

for (name, params, input) in test_cases
    f = () -> mfeval(params, input, M111)
    benchmark_function(f, name, 1, "single_point")
end

## ========================================================================
## BENCHMARK 2: Batch Evaluation (Various Sizes)
## ========================================================================
println("\n--- Benchmark 2: Batch Evaluation (MF 6.1) ---")

batch_sizes = [10, 50, 100, 500, 1000, 2000, 5000, 10000]

for batch_size in batch_sizes
    # Generate batch input: lateral slip sweep
    alpha_sweep = range(-0.3, 0.3, length=batch_size)
    batch_input = Matrix{Float64}(undef, batch_size, 6)
    for i in 1:batch_size
        batch_input[i, :] = [p61.fnomin, 0.0, alpha_sweep[i], 0.0, 0.0, 0.0]
    end
    
    # Pre-allocate output matrix for in-place version
    out_matrix = Matrix{Float64}(undef, batch_size, 30)
    
    f = () -> mfeval!(out_matrix, p61, batch_input, M111)
    
    # Adjust benchmark parameters for larger batches
    evals = batch_size ≤ 1000 ? 100 : (batch_size ≤ 5000 ? 20 : 10)
    samples = batch_size ≤ 1000 ? 500 : (batch_size ≤ 5000 ? 100 : 50)
    
    mean_time = benchmark_function(f, "MF61_batch_N$batch_size", batch_size, "batch"; 
                                  evals=evals, samples=samples)
    
    per_eval = mean_time / batch_size
    @printf("    Per evaluation: %.2f μs\n", per_eval)
end

## ========================================================================
## BENCHMARK 3: Specific Tire Scenarios
## ========================================================================
println("\n--- Benchmark 3: Specific Tire Scenarios ---")

# Helper functions to create test scenarios
function create_combined_slip_grid(Fnom::Float64)
    Nk, Na = 15, 15
    k_vals = range(-0.3, 0.3, length=Nk)
    a_vals = range(-0.3, 0.3, length=Na)
    N = Nk * Na
    grid = Matrix{Float64}(undef, N, 6)
    idx = 1
    for k in k_vals, a in a_vals
        grid[idx, :] = [Fnom, k, a, 0.0, 0.0, 0.0]
        idx += 1
    end
    return grid
end

function create_multi_load_carpet(Fnom::Float64)
    Fz_levels = [0.5, 1.0, 1.5, 2.0] .* Fnom
    N_alpha = 50
    N_total = length(Fz_levels) * N_alpha
    carpet = Matrix{Float64}(undef, N_total, 6)
    alpha_sweep = range(-0.3, 0.3, length=N_alpha)
    
    row = 1
    for Fz in Fz_levels
        for alpha in alpha_sweep
            carpet[row, :] = [Fz, 0.0, alpha, 0.0, 0.0, 0.0]
            row += 1
        end
    end
    return carpet
end

scenarios = [
    ("pure_lateral_sweep", 
     [fill(p61.fnomin, 200) zeros(200) range(-0.5, 0.5, length=200) zeros(200, 3)]),
    ("pure_longitudinal_sweep", 
     [fill(p61.fnomin, 200) range(-0.5, 0.5, length=200) zeros(200, 4)]),
    ("combined_slip_grid", 
     create_combined_slip_grid(p61.fnomin)),
    ("load_sweep", 
     [range(0.5*p61.fnomin, 2.0*p61.fnomin, length=100) zeros(100) fill(0.1, 100) zeros(100, 3)]),
    ("camber_sweep", 
     [fill(p61.fnomin, 50) zeros(50, 2) range(-p61.cammax, p61.cammax, length=50) zeros(50, 2)]),
    ("multi_load_carpet", 
     create_multi_load_carpet(p61.fnomin))
]

for (name, input_matrix) in scenarios
    N = size(input_matrix, 1)
    out_matrix = Matrix{Float64}(undef, N, 30)
    
    f = () -> mfeval!(out_matrix, p61, input_matrix, M111)
    
    # Adjust benchmark parameters based on size
    evals = N ≤ 200 ? 50 : (N ≤ 500 ? 20 : 10)
    samples = N ≤ 200 ? 200 : (N ≤ 500 ? 100 : 50)
    
    mean_time = benchmark_function(f, name, N, "scenario"; evals=evals, samples=samples)
    
    per_eval = mean_time / N
    @printf("    Per evaluation: %.2f μs\n", per_eval)
end

## ========================================================================
## BENCHMARK 4: useMode Performance Impact
## ========================================================================
println("\n--- Benchmark 4: useMode Performance Impact ---")

test_input = MFInputs(p61.fnomin, -0.05, 0.1, 0.02, 0.0, p61.longvl)
useModes = [111, 112, 121, 122, 211, 212, 221, 222]

for useMode_test in useModes
    try
        modes = MFModes(useMode_test)
        
        if modes.is_valid
            f = () -> mfeval(p61, test_input, modes)
            benchmark_function(f, "useMode_$useMode_test", 1, "usemode")
        else
            # Test invalid useMode (should return NaN quickly)
            f = () -> mfeval(p61, test_input, modes)
            benchmark_function(f, "useMode_$(useMode_test)_invalid", 1, "usemode")
        end
    catch e
        @printf("  useMode %-3d: ERROR - %s\n", useMode_test, e)
    end
end

## ========================================================================
## BENCHMARK 5: Memory and Cache Effects
## ========================================================================
println("\n--- Benchmark 5: Cache and Memory Effects ---")

N_cache_test = 1000

# Same input repeated (cache friendly)
same_input = MFInputs(p61.fnomin, -0.1, 0.1, 0.0, 0.0, p61.longvl)
f_same = () -> begin
    for _ in 1:N_cache_test
        mfeval(p61, same_input, M111)
    end
end

b_same = @benchmark ($f_same()) evals=10 samples=100
same_input_time = median(b_same).time / N_cache_test / 1e3  # μs per call

# Different inputs each time (cache unfriendly)
Random.seed!(12345)  # For reproducibility
different_inputs = [MFInputs(p61.fnomin, 
                            rand() * 0.2 - 0.1,  # kappa
                            rand() * 0.2 - 0.1,  # alpha  
                            0.0, 0.0, p61.longvl) for _ in 1:N_cache_test]

f_different = () -> begin
    for inp in different_inputs
        mfeval(p61, inp, M111)
    end
end

b_different = @benchmark ($f_different()) evals=10 samples=100  
different_input_time = median(b_different).time / N_cache_test / 1e3  # μs per call

@printf("  Same input repeated:     %.2f μs per call\n", same_input_time)
@printf("  Different inputs:        %.2f μs per call\n", different_input_time)
@printf("  Cache penalty:           %.1f%%\n", (different_input_time / same_input_time - 1) * 100)

# Store results
push!(results, (test_name="cache_same_input", n_points=1, mean_time_us=same_input_time,
               min_time_us=NaN, max_time_us=NaN, std_time_us=NaN, p95_time_us=NaN, category="cache"))
push!(results, (test_name="cache_different_inputs", n_points=1, mean_time_us=different_input_time,
               min_time_us=NaN, max_time_us=NaN, std_time_us=NaN, p95_time_us=NaN, category="cache"))

## ========================================================================
## BENCHMARK 6: Threading Performance (Julia-specific)
## ========================================================================
println("\n--- Benchmark 6: Threading Performance ---")

if Threads.nthreads() > 1
    # Large batch for threading test
    N_thread = 10000
    thread_input = Matrix{Float64}(undef, N_thread, 6)
    alpha_vals = range(-0.3, 0.3, length=N_thread)
    for i in 1:N_thread
        thread_input[i, :] = [p61.fnomin, 0.0, alpha_vals[i], 0.0, 0.0, 0.0]
    end
    
    out_matrix = Matrix{Float64}(undef, N_thread, 30)
    
    f_threaded = () -> mfeval!(out_matrix, p61, thread_input, M111)
    
    mean_time = benchmark_function(f_threaded, "threading_N$N_thread", N_thread, "threading"; 
                                  evals=5, samples=20)
    
    per_eval = mean_time / N_thread
    @printf("    Per evaluation: %.2f μs (with %d threads)\n", per_eval, Threads.nthreads())
else
    @printf("  Single-threaded execution (no threading benchmark)\n")
end

## ========================================================================
## BENCHMARK 7: Memory Allocation Analysis  
## ========================================================================
println("\n--- Benchmark 7: Memory Allocation Analysis ---")

# Test allocating vs non-allocating versions
test_input = [p61.fnomin 0.0 0.1 0.0 0.0 p61.longvl]
out_matrix = Matrix{Float64}(undef, 1, 30)

# Allocating version
f_alloc = () -> mfeval(p61, test_input, M111)
b_alloc = @benchmark ($f_alloc())
alloc_time = median(b_alloc).time / 1e3
alloc_bytes = b_alloc.memory

# Non-allocating version  
f_nonalloc = () -> mfeval!(out_matrix, p61, test_input, M111)
b_nonalloc = @benchmark ($f_nonalloc())
nonalloc_time = median(b_nonalloc).time / 1e3
nonalloc_bytes = b_nonalloc.memory

@printf("  Allocating version:      %.2f μs, %d bytes allocated\n", alloc_time, alloc_bytes)
@printf("  Non-allocating version:  %.2f μs, %d bytes allocated\n", nonalloc_time, nonalloc_bytes)
@printf("  Allocation overhead:     %.1f%%\n", (alloc_time / nonalloc_time - 1) * 100)

# Store results
push!(results, (test_name="memory_allocating", n_points=1, mean_time_us=alloc_time,
               min_time_us=NaN, max_time_us=NaN, std_time_us=NaN, p95_time_us=NaN, category="memory"))
push!(results, (test_name="memory_nonallocating", n_points=1, mean_time_us=nonalloc_time,
               min_time_us=NaN, max_time_us=NaN, std_time_us=NaN, p95_time_us=NaN, category="memory"))

## ========================================================================
## Write Results to CSV
## ========================================================================
println("\n--- Writing Results ---")

open(RESULTS_FILE, "w") do io
    # Header
    println(io, "test_name,n_points,mean_time_us,min_time_us,max_time_us,std_time_us,p95_time_us,category")
    
    # Data
    for r in results
        @printf(io, "%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%s\n",
                r.test_name, r.n_points, r.mean_time_us, r.min_time_us, r.max_time_us,
                r.std_time_us, r.p95_time_us, r.category)
    end
end

@printf("Results written to: %s\n", RESULTS_FILE)
@printf("Total tests completed: %d\n", length(results))

## ========================================================================
## Summary Statistics
## ========================================================================
println("\n" * "="^65)
println("  Julia Performance Benchmark Summary")
println("="^65)

# Find representative results
function find_result(name::String)
    for r in results
        if r.test_name == name
            return r
        end
    end
    return nothing
end

single_61 = find_result("MF61_single")
batch_1000 = find_result("MF61_batch_N1000")
lateral_sweep = find_result("pure_lateral_sweep")

if single_61 !== nothing
    @printf("Single point (MF 6.1):       %.2f μs\n", single_61.mean_time_us)
end
if batch_1000 !== nothing
    per_eval = batch_1000.mean_time_us / batch_1000.n_points
    @printf("Batch N=1000 (per eval):     %.2f μs\n", per_eval)
end
if lateral_sweep !== nothing
    per_eval = lateral_sweep.mean_time_us / lateral_sweep.n_points
    @printf("Lateral sweep (N=200):       %.2f μs per point\n", per_eval)
end

println("\nBenchmark completed successfully!")
println("="^65)