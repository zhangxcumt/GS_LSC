%% RUN_TABLE4_SWEEP_COST  Cost of one sweep with and without the cached matrix.
%
%   Table 4 of the paper compares two ways of advancing the block Gauss-Seidel
%   iteration by one sweep:
%
%     (a) with the normal-equation matrix A cached: one matrix-vector product,
%     (b) with the covariance recomputed: N^2 covariance elements evaluated
%         at every sweep.
%
%   (a) is timed in RUN_STANDARD_LSC, which stores the time of one matvec with
%   the assembled matrix in R.info.tMatvec.  (b) is timed here by running the
%   block Gauss-Seidel solver for a single sweep on the same data.
%
%   Nothing is factorised, so this script needs no large matrix; but it does
%   evaluate N^2 covariances, which takes about 30 s at N = 5 x 10^4.
%
%   See also RUN_STANDARD_LSC, LSC_MATVEC.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Table 4: cost of one sweep ===\n');

E = setup_experiment(cfg);

% ---- (a) one matvec, timed during the standard LSC run -------------------
tMatvec = NaN;
fname = fullfile(cfg.resDir, 'standard_lsc.mat');
if exist(fname, 'file')
    S = load(fname, 'R');
    if isfield(S.R.info, 'tMatvec'), tMatvec = S.R.info.tMatvec; end
end
if isnan(tMatvec)
    fprintf('   (a) not available: run run_standard_lsc.m first (it times one matvec).\n');
else
    fprintf('   (a) one sweep with A cached       : %7.3f s\n', tMatvec);
end

% ---- (b) one sweep recomputing the covariance ----------------------------
fprintf('   (b) timing one sweep with the covariance recomputed ...\n');
[~, ~, iOne] = lsc_block_gauss_seidel(E.y, E.ctx, E.dvar, struct( ...
    'blockSize', cfg.blockSize, 'stopRule', 'sweeps', 'nSweep', 1, ...
    'verbose', false));
tSweep = iOne.meanSweepTime;
fprintf('   (b) one sweep with the covariance : %7.3f s\n', tSweep);

if ~isnan(tMatvec)
    fprintf('\n   ratio (b)/(a) = %.1f\n', tSweep / tMatvec);
end
fprintf('   N = %d, covariance elements per sweep = %.3e\n', E.ctx.N, E.ctx.N^2);
