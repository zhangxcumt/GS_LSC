%% RUN_BLOCK_GS_LSC  Block Gauss-Seidel LSC of the sub-region (Table 3).
%
%   Runs the iterative solver for the block sizes of Table 3 (nb = 2000 and
%   4000) and writes results/block_gs_lsc_<nb>.mat for each of them.  The
%   iteration stops on the predicted quantity: when the largest change of the
%   height anomaly at the 222 check points between two sweeps is below 0.01 cm.
%
%   Output (results/block_gs_lsc_<nb>.mat)
%     R.method, R.nb, R.alpha, R.info, R.y, R.ym, R.sel, R.dvar, R.stats,
%     R.selfStd, R.dThetaHist (for Fig. 5)
%
%   See also RUN_STANDARD_LSC, RUN_SEQUENTIAL_LSC, MAKE_TABLE3, MAKE_FIG5.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Block Gauss-Seidel LSC, Colorado sub-region ===\n');

E = setup_experiment(cfg);

% stopping rule of the paper: change of the predicted height anomaly < 0.01 cm
xyzG = BLh_XYZR(E.D.gnssBLh);
thetaFun = @(a) 100 * signal_predict(E.ctx.Gcross, a, E.ctx, ...
    xyzG(:, 1:3), xyzG(:, 4), struct('gnorm', E.D.gnssGnorm));

for nb = [2000 4000]

    fprintf('\n   ---- nb = %d ----\n', nb);
    opt = struct('blockSize', nb, 'order', 'lat', 'stopRule', 'pred', ...
                 'tolPred', 0.01, 'thetaFun', thetaFun, 'nSweep', 300, ...
                 'single', false, 'verbose', true);

    [aBgs, ~, iBgs] = lsc_block_gauss_seidel(E.y, E.ctx, E.dvar, opt);

    [~, S] = predict_gnss(E.ctx.Gcross, aBgs, E.ctx, E.D);
    selfStd = std(-E.dvar .* aBgs);

    R = struct('method', 'blockgs', 'nb', nb, 'alpha', aBgs, 'info', iBgs, ...
               'y', E.y, 'ym', E.ym, 'sel', E.sel, 'dvar', E.dvar, ...
               'n', E.ctx.N, 'stats', S, 'selfStd', selfStd, 'cov', cfg.cov, ...
               'dThetaHist', iBgs.dThetaHist);

    fname = fullfile(cfg.resDir, sprintf('block_gs_lsc_%d.mat', nb));
    save(fname, 'R', '-v7.3');
    fprintf('   saved %s\n', fname);

    fprintf('\n   nb = %d   N = %d   RMS = %.4f cm   STD = %.4f cm\n', ...
        nb, R.n, S.rms * 100, S.std * 100);
    fprintf('   sweeps = %d (%s)   time = %.1f s (%.2f s per sweep)\n', ...
        iBgs.sweeps, iBgs.stopReason, iBgs.time, iBgs.meanSweepTime);
    fprintf('   residual at stop = %.3e   storage = %.0f MB\n', ...
        iBgs.relResid, iBgs.blockBytes / 2^20);
end
