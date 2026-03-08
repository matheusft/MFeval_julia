% performance_benchmark_matlab.m
%
% Comprehensive performance benchmarking suite for MATLAB mfeval
% Generates timing data across various scenarios for comparison with Julia
%
% Prerequisites:
%   1. MATLAB with mfeval toolbox on path
%   2. TIR files in test/fixtures/
%
% Output:
%   test/matlab_performance_results.csv - Detailed timing data
%
% Usage:
%   cd /path/to/MFeval_julia
%   matlab -batch "run('test/performance_benchmark_matlab.m')"

fprintf('=================================================================\n');
fprintf('  MATLAB mfeval Performance Benchmark Suite\n');
fprintf('  Platform: MATLAB %s\n', version);
fprintf('  Date: %s\n', datestr(now));
fprintf('=================================================================\n\n');

% Configuration
FIX_DIR = fullfile('test', 'fixtures');
RESULTS_FILE = fullfile('test', 'matlab_performance_results.csv');

% Check if fixtures exist
if ~exist(FIX_DIR, 'dir')
    error('Fixtures directory not found: %s', FIX_DIR);
end

% Load TIR files
fprintf('Loading TIR files...\n');
tir61_path = fullfile(FIX_DIR, 'MagicFormula61_Parameters.tir');
tir62_path = fullfile(FIX_DIR, 'MagicFormula62_Parameters.tir');
tir52_path = fullfile(FIX_DIR, 'MagicFormula52_Parameters.tir');

if ~exist(tir61_path, 'file'), error('MF61 TIR file not found'); end
if ~exist(tir62_path, 'file'), error('MF62 TIR file not found'); end  
if ~exist(tir52_path, 'file'), error('MF52 TIR file not found'); end

p61 = mfeval.readTIR(tir61_path);
p62 = mfeval.readTIR(tir62_path);
p52 = mfeval.readTIR(tir52_path);

fprintf('  MF 6.1: Fnom=%.1f N, Vnom=%.1f m/s\n', p61.FNOMIN, p61.LONGVL);
fprintf('  MF 6.2: Fnom=%.1f N, Vnom=%.1f m/s\n', p62.FNOMIN, p62.LONGVL);
fprintf('  MF 5.2: Fnom=%.1f N, Vnom=%.1f m/s\n', p52.FNOMIN, p52.LONGVL);

% Benchmark configuration
useMode = 111;
N_WARMUP = 100;      % Warmup iterations
N_TIMING = 1000;     % Timing iterations for single point tests
N_BATCH_TIMING = 10; % Timing iterations for batch tests

% Results storage
results = {};
result_idx = 1;

%% ========================================================================
%% BENCHMARK 1: Single Point Evaluation (All MF Versions)
%% ========================================================================
fprintf('\n--- Benchmark 1: Single Point Evaluation ---\n');

test_cases = {
    {'MF61_single', tir61_path, [p61.FNOMIN, -0.05, 0.1, 0.02, 0.0, p61.LONGVL]},
    {'MF62_single', tir62_path, [p62.FNOMIN, -0.05, 0.1, 0.02, 0.0, p62.LONGVL]},
    {'MF52_single', tir52_path, [p52.FNOMIN, -0.05, 0.1, 0.02, 0.0, p52.LONGVL]}
};

for i = 1:length(test_cases)
    name = test_cases{i}{1};
    tir_path = test_cases{i}{2};
    input = test_cases{i}{3};
    
    % Warmup
    for w = 1:N_WARMUP
        mfeval(tir_path, input, useMode);
    end
    
    % Timing
    times = zeros(N_TIMING, 1);
    for t = 1:N_TIMING
        tic;
        mfeval(tir_path, input, useMode);
        times(t) = toc;
    end
    
    % Statistics
    mean_time = mean(times) * 1e6;  % Convert to microseconds
    min_time = min(times) * 1e6;
    max_time = max(times) * 1e6;
    std_time = std(times) * 1e6;
    p95_time = prctile(times, 95) * 1e6;
    
    fprintf('  %-15s: %.2f ± %.2f μs (min=%.2f, max=%.2f, p95=%.2f)\n', ...
            name, mean_time, std_time, min_time, max_time, p95_time);
    
    % Store result
    results{result_idx} = {name, 1, mean_time, min_time, max_time, std_time, p95_time, 'single_point'};
    result_idx = result_idx + 1;
end

%% ========================================================================
%% BENCHMARK 2: Batch Evaluation (Various Sizes)
%% ========================================================================
fprintf('\n--- Benchmark 2: Batch Evaluation (MF 6.1) ---\n');

batch_sizes = [10, 50, 100, 500, 1000, 2000, 5000, 10000];

for batch_size = batch_sizes
    % Generate batch input: lateral slip sweep
    alpha_sweep = linspace(-0.3, 0.3, batch_size)';
    batch_input = [repmat(p61.FNOMIN, batch_size, 1), ...
                   zeros(batch_size, 1), ...
                   alpha_sweep, ...
                   zeros(batch_size, 3)];
    
    % Warmup
    for w = 1:5
        mfeval(tir61_path, batch_input, useMode);
    end
    
    % Timing
    times = zeros(N_BATCH_TIMING, 1);
    for t = 1:N_BATCH_TIMING
        tic;
        mfeval(tir61_path, batch_input, useMode);
        times(t) = toc;
    end
    
    % Statistics
    mean_time = mean(times) * 1e3;  % Convert to milliseconds
    min_time = min(times) * 1e3;
    max_time = max(times) * 1e3;
    std_time = std(times) * 1e3;
    p95_time = prctile(times, 95) * 1e3;
    per_eval = mean_time * 1e3 / batch_size;  % μs per evaluation
    
    fprintf('  N=%-6d: %.3f ± %.3f ms (%.2f μs/eval, min=%.3f, max=%.3f)\n', ...
            batch_size, mean_time, std_time, per_eval, min_time, max_time);
    
    % Store result
    test_name = sprintf('MF61_batch_N%d', batch_size);
    results{result_idx} = {test_name, batch_size, per_eval, min_time*1e3/batch_size, ...
                          max_time*1e3/batch_size, std_time*1e3/batch_size, ...
                          p95_time*1e3/batch_size, 'batch'};
    result_idx = result_idx + 1;
end

%% ========================================================================
%% BENCHMARK 3: Specific Tire Scenarios
%% ========================================================================
fprintf('\n--- Benchmark 3: Specific Tire Scenarios ---\n');

scenarios = {
    {'pure_lateral_sweep', [repmat(p61.FNOMIN, 200, 1), zeros(200, 1), ...
                           linspace(-0.5, 0.5, 200)', zeros(200, 3)]},
    {'pure_longitudinal_sweep', [repmat(p61.FNOMIN, 200, 1), ...
                                linspace(-0.5, 0.5, 200)', zeros(200, 4)]},
    {'combined_slip_grid', create_combined_slip_grid(p61.FNOMIN)},
    {'load_sweep', [linspace(0.5*p61.FNOMIN, 2.0*p61.FNOMIN, 100)', ...
                   zeros(100, 1), repmat(0.1, 100, 1), zeros(100, 3)]},
    {'camber_sweep', [repmat(p61.FNOMIN, 50, 1), zeros(50, 2), ...
                     linspace(-p61.CAMMAX, p61.CAMMAX, 50)', zeros(50, 2)]},
    {'multi_load_carpet', create_multi_load_carpet(p61.FNOMIN)}
};

for s = 1:length(scenarios)
    name = scenarios{s}{1};
    input_matrix = scenarios{s}{2};
    N = size(input_matrix, 1);
    
    % Warmup
    for w = 1:3
        mfeval(tir61_path, input_matrix, useMode);
    end
    
    % Timing
    times = zeros(N_BATCH_TIMING, 1);
    for t = 1:N_BATCH_TIMING
        tic;
        mfeval(tir61_path, input_matrix, useMode);
        times(t) = toc;
    end
    
    % Statistics
    mean_time = mean(times) * 1e3;  % ms
    min_time = min(times) * 1e3;
    max_time = max(times) * 1e3;
    std_time = std(times) * 1e3;
    p95_time = prctile(times, 95) * 1e3;
    per_eval = mean_time * 1e3 / N;  % μs per evaluation
    
    fprintf('  %-25s (N=%4d): %.3f ms (%.2f μs/eval)\n', ...
            name, N, mean_time, per_eval);
    
    % Store result
    results{result_idx} = {name, N, per_eval, min_time*1e3/N, ...
                          max_time*1e3/N, std_time*1e3/N, ...
                          p95_time*1e3/N, 'scenario'};
    result_idx = result_idx + 1;
end

%% ========================================================================
%% BENCHMARK 4: useMode Performance Impact
%% ========================================================================
fprintf('\n--- Benchmark 4: useMode Performance Impact ---\n');

test_input = [p61.FNOMIN, -0.05, 0.1, 0.02, 0.0, p61.LONGVL];
useModes = [111, 112, 121, 122, 211, 212, 221, 222];

for useMode_test = useModes
    % Warmup
    for w = 1:N_WARMUP
        try
            mfeval(tir61_path, test_input, useMode_test);
        catch
            % Skip invalid useModes
            continue;
        end
    end
    
    % Timing
    times = zeros(N_TIMING, 1);
    valid_times = 0;
    for t = 1:N_TIMING
        try
            tic;
            mfeval(tir61_path, test_input, useMode_test);
            times(valid_times + 1) = toc;
            valid_times = valid_times + 1;
        catch
            % Skip invalid useModes
            continue;
        end
    end
    
    if valid_times == 0
        fprintf('  useMode %-3d: INVALID\n', useMode_test);
        continue;
    end
    
    % Statistics
    times = times(1:valid_times);
    mean_time = mean(times) * 1e6;  % μs
    std_time = std(times) * 1e6;
    
    fprintf('  useMode %-3d: %.2f ± %.2f μs\n', useMode_test, mean_time, std_time);
    
    % Store result
    test_name = sprintf('useMode_%d', useMode_test);
    results{result_idx} = {test_name, 1, mean_time, min(times)*1e6, ...
                          max(times)*1e6, std_time, prctile(times, 95)*1e6, 'usemode'};
    result_idx = result_idx + 1;
end

%% ========================================================================
%% BENCHMARK 5: Memory and Cache Effects  
%% ========================================================================
fprintf('\n--- Benchmark 5: Cache and Memory Effects ---\n');

% Test effect of repeated evaluations with same vs different inputs
N_cache_test = 1000;

% Same input repeated (cache friendly)
same_input = [p61.FNOMIN, -0.1, 0.1, 0.0, 0.0, p61.LONGVL];
tic;
for i = 1:N_cache_test
    mfeval(tir61_path, same_input, useMode);
end
same_input_time = toc / N_cache_test * 1e6;  % μs per call

% Different inputs each time (cache unfriendly)  
different_inputs = [repmat(p61.FNOMIN, N_cache_test, 1), ...
                   rand(N_cache_test, 1) * 0.2 - 0.1, ...  % kappa
                   rand(N_cache_test, 1) * 0.2 - 0.1, ...  % alpha
                   zeros(N_cache_test, 3)];

tic;
for i = 1:N_cache_test
    mfeval(tir61_path, different_inputs(i, :), useMode);
end
different_input_time = toc / N_cache_test * 1e6;  % μs per call

fprintf('  Same input repeated:     %.2f μs per call\n', same_input_time);
fprintf('  Different inputs:        %.2f μs per call\n', different_input_time);
fprintf('  Cache penalty:           %.1f%%\n', (different_input_time / same_input_time - 1) * 100);

% Store results
results{result_idx} = {'cache_same_input', 1, same_input_time, NaN, NaN, NaN, NaN, 'cache'};
result_idx = result_idx + 1;
results{result_idx} = {'cache_different_inputs', 1, different_input_time, NaN, NaN, NaN, NaN, 'cache'};
result_idx = result_idx + 1;

%% ========================================================================
%% Write Results to CSV
%% ========================================================================
fprintf('\n--- Writing Results ---\n');

fid = fopen(RESULTS_FILE, 'w');
if fid == -1
    error('Could not open results file: %s', RESULTS_FILE);
end

% Header
fprintf(fid, 'test_name,n_points,mean_time_us,min_time_us,max_time_us,std_time_us,p95_time_us,category\n');

% Data
for i = 1:length(results)
    r = results{i};
    fprintf(fid, '%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%s\n', ...
            r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, r{8});
end

fclose(fid);

fprintf('Results written to: %s\n', RESULTS_FILE);
fprintf('Total tests completed: %d\n', length(results));

%% ========================================================================
%% Summary Statistics
%% ========================================================================
fprintf('\n=================================================================\n');
fprintf('  MATLAB Performance Benchmark Summary\n');
fprintf('=================================================================\n');

% Find representative results
single_61 = find_result(results, 'MF61_single');
batch_1000 = find_result(results, 'MF61_batch_N1000');
lateral_sweep = find_result(results, 'pure_lateral_sweep');

if ~isempty(single_61)
    fprintf('Single point (MF 6.1):       %.2f μs\n', single_61{3});
end
if ~isempty(batch_1000)  
    fprintf('Batch N=1000 (per eval):     %.2f μs\n', batch_1000{3});
end
if ~isempty(lateral_sweep)
    fprintf('Lateral sweep (N=200):       %.2f μs per point\n', lateral_sweep{3});
end

fprintf('\nBenchmark completed successfully!\n');
fprintf('=================================================================\n');

%% Helper Functions
function grid = create_combined_slip_grid(Fnom)
    Nk = 15; Na = 15;
    [KK, AA] = meshgrid(linspace(-0.3, 0.3, Nk), linspace(-0.3, 0.3, Na));
    kv = KK(:); av = AA(:);
    N = length(kv);
    grid = [repmat(Fnom, N, 1), kv, av, zeros(N, 3)];
end

function carpet = create_multi_load_carpet(Fnom)
    Fz_levels = [0.5, 1.0, 1.5, 2.0] * Fnom;
    N_alpha = 50;
    N_total = length(Fz_levels) * N_alpha;
    carpet = zeros(N_total, 6);
    row = 1;
    for Fz = Fz_levels
        alpha_sweep = linspace(-0.3, 0.3, N_alpha)';
        for i = 1:N_alpha
            carpet(row, :) = [Fz, 0, alpha_sweep(i), 0, 0, 0];
            row = row + 1;
        end
    end
end

function result = find_result(results, name)
    result = {};
    for i = 1:length(results)
        if strcmp(results{i}{1}, name)
            result = results{i};
            break;
        end
    end
end