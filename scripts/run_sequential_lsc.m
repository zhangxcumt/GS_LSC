%% RUN_SEQUENTIAL_LSC  Sequential (stepwise) LSC of the sub-region (Table 3).
%
%   Moritz/Tscherning recursion: only n_b x n_b matrices are inverted, but the
%   recursion state P is N x N in double precision (20 GB at N = 50 902).  The
%   function refuses to run when P exceeds opt.maxMemGB.
%
%   Output (results/sequential_lsc.mat)
%     R.method, R.alpha, R.info, R.y, R.ym, R.sel, R.dvar, R.stats, R.selfStd
%
%   See also RUN_STANDARD_LSC, RUN_BLOCK_GS_LSC, MAKE_TABLE3.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Sequential LSC, Colorado sub-region ===\n');

% Precision of the recursion state.  P is N x N, so it needs 8 N^2 bytes in
% double precision (19.3 GiB at N = 50 902, the configuration of the paper) and
% 4 N^2 bytes in single precision.  Single precision costs a relative accuracy
% of about 1e-6 in the coefficients, far below what a height anomaly in cm
% needs, and it is the only way to run the full sub-region on a machine where
% the double-precision state does not fit.  Set the environment variable
% GSLSC_SEQ_SINGLE=1 to use it.
singlePrecision = strcmpi(getenv('GSLSC_SEQ_SINGLE'), '1');
if singlePrecision
    fprintf('   [single precision recursion state, set by GSLSC_SEQ_SINGLE]\n');
end

E = setup_experiment(cfg);

opt = struct('blockSize', cfg.blockSize, 'order', 'lat', 'single', singlePrecision, ...
             'maxMemGB', 26, 'verbose', true);

[aSeq, ~, iSeq] = lsc_sequential(E.y, E.ctx, E.dvar, opt);

[~, S] = predict_gnss(E.ctx.Gcross, aSeq, E.ctx, E.D);
selfStd = std(-E.dvar .* aSeq);

R = struct('method', 'sequential', 'alpha', aSeq, 'info', iSeq, ...
           'y', E.y, 'ym', E.ym, 'sel', E.sel, 'dvar', E.dvar, ...
           'n', E.ctx.N, 'stats', S, 'selfStd', selfStd, 'cov', cfg.cov);

fname = fullfile(cfg.resDir, 'sequential_lsc.mat');
save(fname, 'R', '-v7.3');
fprintf('   saved %s\n', fname);

fprintf('\n   N = %d   RMS = %.4f cm   STD = %.4f cm   self-fit STD = %.4f mGal\n', ...
    R.n, S.rms * 100, S.std * 100, selfStd);
fprintf('   time = %.1f s   recursion state P = %.2f GB\n', iSeq.time, iSeq.pBytes / 2^30);
