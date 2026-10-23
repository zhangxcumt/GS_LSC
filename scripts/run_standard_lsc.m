%% RUN_STANDARD_LSC  Standard LSC of the paper's sub-region (Table 3).
%
%   Solves the collocation equations with the Cholesky factorisation of the
%   full N x N matrix A = C + D, predicts the residual height anomaly at the
%   222 GNSS/levelling check points, and writes results/standard_lsc.mat.
%
%   At N = 50 902 the matrix needs 9.65 GB in single precision and the
%   factorisation about 20 GB of peak memory.  Set the environment variable
%   GSLSC_NMAX (see SETUP_EXPERIMENT) to run a reduced experiment first.
%
%   Output (results/standard_lsc.mat)
%     R.method, R.alpha, R.info, R.y, R.ym, R.sel, R.dvar, R.stats, R.selfStd
%
%   See also RUN_SEQUENTIAL_LSC, RUN_BLOCK_GS_LSC, MAKE_TABLE3.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Standard LSC, Colorado sub-region ===\n');

E = setup_experiment(cfg);

opt = struct('precision', cfg.precision, 'blockSize', cfg.blockSize, ...
             'verbose', true, 'timeMatvec', true);

[aStd, iStd] = lsc_standard(E.y, E.ctx, E.dvar, opt);

[~, S] = predict_gnss(E.ctx.Gcross, aStd, E.ctx, E.D);
selfStd = std(-E.dvar .* aStd);

R = struct('method', 'standard', 'alpha', aStd, 'info', iStd, ...
           'y', E.y, 'ym', E.ym, 'sel', E.sel, 'dvar', E.dvar, ...
           'n', E.ctx.N, 'stats', S, 'selfStd', selfStd, 'cov', cfg.cov);

fname = fullfile(cfg.resDir, 'standard_lsc.mat');
save(fname, 'R', '-v7.3');
fprintf('   saved %s\n', fname);

fprintf('\n   N = %d   RMS = %.4f cm   STD = %.4f cm   self-fit STD = %.4f mGal\n', ...
    R.n, S.rms * 100, S.std * 100, selfStd);
fprintf('   time = %.1f s   storage = %.2f GB\n', iStd.time, iStd.bytes / 2^30);
