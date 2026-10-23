%% MAKE_TABLE3  Assemble Table 3 of the paper from the saved solutions.
%
%   Reads the result files written by
%
%       run_standard_lsc.m      results/standard_lsc.mat
%       run_sequential_lsc.m    results/sequential_lsc.mat
%       run_block_gs_lsc.m      results/block_gs_lsc_<nb>.mat
%
%   and prints the modelling accuracy, the run time and the storage of each
%   solver, in the layout of Table 3.  Missing files are reported as "-", so
%   the table can be built incrementally.
%
%   The storage figures are the analytic requirements of the three algorithms
%   (see STORAGE_MODEL): 4 N^2 bytes for the single-precision normal-equation
%   matrix of standard LSC, 8 N^2 bytes for the double-precision recursion
%   state of sequential LSC, and 8 nb^2 bytes for one diagonal block of the
%   block Gauss-Seidel solver.
%
%   A CSV copy is written to results/table3.csv.
%
%   See also STORAGE_MODEL, MAKE_FIG1_STORAGE.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Table 3: modelling accuracy, time and storage ===\n');

files = { 'standard_lsc.mat',   'Standard LSC',                   'standard',   [];
          'sequential_lsc.mat', 'Sequential LSC',                 'sequential', [];
          'block_gs_lsc_2000.mat', 'Block Gauss-Seidel (nb = 2000)', 'blockgs', 2000;
          'block_gs_lsc_4000.mat', 'Block Gauss-Seidel (nb = 4000)', 'blockgs', 4000 };

fprintf('\n   %-32s %9s %9s %9s %9s %12s\n', 'method', 'RMS(cm)', 'STD(cm)', ...
    'selfSTD', 'time(s)', 'storage');
fprintf('   %s\n', repmat('-', 1, 86));

rows = {};
for k = 1:size(files, 1)
    fname = fullfile(cfg.resDir, files{k, 1});
    if ~exist(fname, 'file')
        fprintf('   %-32s %9s %9s %9s %9s %12s\n', files{k, 2}, '-', '-', '-', '-', '-');
        continue;
    end
    S = load(fname, 'R');
    R = S.R;
    n  = R.n;
    nb = files{k, 4};
    if isempty(nb)
        if isfield(R.info, 'blockSize')
            nb = R.info.blockSize;
        else
            nb = cfg.blockSize;
        end
    end

    bytes = storage_model(files{k, 3}, n, nb);
    if strcmp(files{k, 3}, 'standard') && isfield(R.info, 'bytes')
        bytes = R.info.bytes;                   % measured, single precision
    elseif strcmp(files{k, 3}, 'sequential') && isfield(R.info, 'pBytes')
        bytes = R.info.pBytes;                  % measured
    end

    if bytes >= 2^30
        sTxt = sprintf('%.2f GB', bytes / 2^30);
    else
        sTxt = sprintf('%.0f MB', bytes / 2^20);
    end

    name = files{k, 2};
    if isfield(R.info, 'single') && R.info.single
        name = [name ' (P in single)'];     %#ok<AGROW>
    end

    fprintf('   %-32s %9.4f %9.4f %9.4f %9.1f %12s\n', ...
        name, R.stats.rms * 100, R.stats.std * 100, R.selfStd, ...
        R.info.time, sTxt);
    rows(end + 1, :) = {name, R.stats.rms * 100, R.stats.std * 100, ...
                        R.selfStd, R.info.time, sTxt};               %#ok<SAGROW>
end

if ~isempty(rows)
    T = cell2table(rows, 'VariableNames', ...
        {'method', 'rms_cm', 'std_cm', 'selffit_std_mGal', 'time_s', 'storage'});
    csv = fullfile(cfg.resDir, 'table3.csv');
    writetable(T, csv);
    fprintf('\n   written %s\n', csv);
end

fprintf(['\n   Storage is reported for the actual N of the run; the paper quotes\n' ...
         '   9.31 GB, 18.63 GB and 32 MB / 128 MB for N = 5 x 10^4.\n']);
