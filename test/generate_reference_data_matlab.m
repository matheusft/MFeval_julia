% generate_reference_data_matlab.m
%
% Generates the CSV reference files for MFeval.jl Phase 4 regression tests
% using the original MATLAB mfeval toolbox as the external ground truth.
%
% Prerequisites
% -------------
%   1. MATLAB with the mfeval toolbox on the path:
%        addpath('/path/to/mfeval')
%   2. Run from the MFeval.jl/ root directory (or adjust FIX_DIR / OUT_DIR).
%
% Output
% ------
%   test/reference_data/*.csv  — 30-column CSV files, header row included.
%   Each file replaces the Julia-generated version, enabling the Phase 4
%   regression tests to validate the Julia solver against true MATLAB output.
%
% Usage
% -----
%   cd /path/to/MFeval.jl
%   matlab -batch "run('test/generate_reference_data_matlab.m')"

FIX_DIR = fullfile('test', 'fixtures');
OUT_DIR = fullfile('test', 'reference_data');

if ~exist(OUT_DIR, 'dir'); mkdir(OUT_DIR); end

COLS = {'Fx','Fy','Fz','Mx','My','Mz', ...
        'kappa','alpha','gamma','phit','Vx','pressure', ...
        'Re','rho','two_a','t','mux','muy','omega','Rl','two_b','Mzr', ...
        'Cx','Cy','Cz','Kya','sigmax','sigmay','inst_Kya','Kxk'};

tir61 = fullfile(FIX_DIR, 'MagicFormula61_Parameters.tir');
tir62 = fullfile(FIX_DIR, 'MagicFormula62_Parameters.tir');
tir52 = fullfile(FIX_DIR, 'MagicFormula52_Parameters.tir');

p61 = mfeval.readTIR(tir61);
p62 = mfeval.readTIR(tir62);
p52 = mfeval.readTIR(tir52);

Fnom61 = p61.FNOMIN;   Vnom61 = p61.LONGVL;
Fnom62 = p62.FNOMIN;   Vnom62 = p62.LONGVL;
Fnom52 = p52.FNOMIN;   Vnom52 = p52.LONGVL;

useMode = 111;

% ── MF 6.1 ────────────────────────────────────────────────────────────────────

N = 200;
alpha = linspace(-0.5, 0.5, N)';
inp = [repmat(Fnom61,N,1), zeros(N,1), alpha, zeros(N,3)];
write_csv(fullfile(OUT_DIR,'mf61_pure_lat.csv'), mfeval(tir61,inp,useMode), COLS);

kappa = linspace(-0.5, 0.5, N)';
inp = [repmat(Fnom61,N,1), kappa, zeros(N,4)];
write_csv(fullfile(OUT_DIR,'mf61_pure_lon.csv'), mfeval(tir61,inp,useMode), COLS);

Nk = 11; Na = 11;
[KK,AA] = meshgrid(linspace(-0.3,0.3,Nk), linspace(-0.3,0.3,Na));
kv = KK(:); av = AA(:); N2 = numel(kv);
inp = [repmat(Fnom61,N2,1), kv, av, zeros(N2,3)];
write_csv(fullfile(OUT_DIR,'mf61_combined.csv'), mfeval(tir61,inp,useMode), COLS);

N = 50;
gamma = linspace(-p61.CAMMAX, p61.CAMMAX, N)';
inp = [repmat(Fnom61,N,1), zeros(N,1), zeros(N,1), gamma, zeros(N,2)];
write_csv(fullfile(OUT_DIR,'mf61_camber.csv'), mfeval(tir61,inp,useMode), COLS);

Fz = linspace(0.25*Fnom61, 2.0*Fnom61, N)';
inp = [Fz, zeros(N,1), repmat(0.1,N,1), zeros(N,3)];
write_csv(fullfile(OUT_DIR,'mf61_load.csv'), mfeval(tir61,inp,useMode), COLS);

N = 30;
pres = linspace(p61.PRESMIN*1.1, p61.PRESMAX*0.9, N)';
inp = [repmat(Fnom61,N,1), zeros(N,1), repmat(0.1,N,1), zeros(N,2), repmat(Vnom61,N,1), pres];
write_csv(fullfile(OUT_DIR,'mf61_pressure.csv'), mfeval(tir61,inp,useMode), COLS);

% Edge cases
edge_inp = [
    0,          0,   0,           0,   0,  Vnom61;   % Fz=0
    Fnom61,     0,   0,           0,   0,  0;         % standstill
    Fnom61,    -1,   0,           0,   0,  Vnom61;   % locked wheel
    Fnom61,     0,   0,  p61.CAMMAX,   0,  Vnom61;   % +max camber
    Fnom61,     0,   0, -p61.CAMMAX,   0,  Vnom61;   % -max camber
    Fnom61,     0,   p61.ALPMAX*1.5, 0, 0, Vnom61;   % alpha > ALPMAX
    Fnom61, p61.KPUMAX*1.5, 0,   0,   0,  Vnom61;   % kappa > KPUMAX
    Fnom61,     0,   0.1,         0,   0,  0.3;       % low speed
    Fnom61*2,   0,   0.1,         0,   0,  Vnom61;   % high load
];
write_csv(fullfile(OUT_DIR,'mf61_edge_cases.csv'), mfeval(tir61,edge_inp,useMode), COLS);

% ── MF 6.2 ────────────────────────────────────────────────────────────────────
N = 200;
alpha = linspace(-0.5, 0.5, N)';
inp = [repmat(Fnom62,N,1), zeros(N,1), alpha, zeros(N,3)];
write_csv(fullfile(OUT_DIR,'mf62_pure_lat.csv'), mfeval(tir62,inp,useMode), COLS);

% ── MF 5.2 ────────────────────────────────────────────────────────────────────
alpha = linspace(-0.5, 0.5, N)';
inp = [repmat(Fnom52,N,1), zeros(N,1), alpha, zeros(N,3)];
write_csv(fullfile(OUT_DIR,'mf52_pure_lat.csv'), mfeval(tir52,inp,useMode), COLS);

fprintf('Done. All reference CSVs written to %s\n', OUT_DIR);

% ── Helper ────────────────────────────────────────────────────────────────────
function write_csv(path, mat, cols)
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', strjoin(cols, ','));
    for i = 1:size(mat,1)
        fprintf(fid, '%.17g', mat(i,1));
        fprintf(fid, ',%.17g', mat(i,2:end));
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('  wrote %s  (%d rows)\n', path, size(mat,1));
end