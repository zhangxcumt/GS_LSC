%% MAKE_FIG5_CONVERGENCE  Fig. 5: convergence of the predicted quantity.
%
%   Reproduces Fig. 5 of the paper: for n_b = 2000 and n_b = 4000, the largest
%   change of the height anomaly at the 222 GNSS/levelling check points between
%   two successive sweeps, with the 0.01 cm stopping threshold.
%
%   The histories are read from results/block_gs_lsc_<nb>.mat, which
%   RUN_BLOCK_GS_LSC writes.  If such a file is missing, the iteration is
%   rerun here, which costs one sweep per point of the curve (about 30 s per
%   sweep at N = 50 902).
%
%   Writes figures/fig5_convergence.png and .pdf, and results/fig5_data.mat.
%
%   See also PLOT_FIG5, LSC_BLOCK_GAUSS_SEIDEL, RUN_BLOCK_GS_LSC.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Fig. 5: convergence of the predicted quantity ===\n');

nbList = [2000 4000];
store  = cell(numel(nbList), 1);
needRun = false;

for k = 1:numel(nbList)
    nb    = nbList(k);
    fname = fullfile(cfg.resDir, sprintf('block_gs_lsc_%d.mat', nb));
    if exist(fname, 'file')
        S = load(fname, 'R');
        hist = S.R.info.dThetaHist;
        store{k} = struct('nb', nb, 'dTheta', hist(:));
        fprintf('   nb = %d: %d sweeps read from %s\n', nb, numel(hist), ...
            sprintf('block_gs_lsc_%d.mat', nb));
    else
        fprintf('   nb = %d: no saved history, the iteration will be rerun\n', nb);
        store{k} = struct('nb', nb, 'dTheta', []);
        needRun  = true;
    end
end

% ---- rerun only what is missing -----------------------------------------
if needRun
    E = setup_experiment(cfg);
    xyzG = BLh_XYZR(E.D.gnssBLh);
    thetaFun = @(a) 100 * signal_predict(E.ctx.Gcross, a, E.ctx, ...
        xyzG(:, 1:3), xyzG(:, 4), struct('gnorm', E.D.gnssGnorm));

    for k = 1:numel(nbList)
        if ~isempty(store{k}.dTheta), continue; end
        fprintf('\n   ---- nb = %d ----\n', nbList(k));
        [~, ~, iBgs] = lsc_block_gauss_seidel(E.y, E.ctx, E.dvar, struct( ...
            'blockSize', nbList(k), 'stopRule', 'sweeps', 'nSweep', 60, ...
            'thetaFun', thetaFun, 'verbose', true));
        store{k}.dTheta = iBgs.dThetaHist;
    end
end

save(fullfile(cfg.resDir, 'fig5_data.mat'), 'store', 'nbList');
plot_fig5(store, nbList, cfg);
fprintf('\n   written %s\n', fullfile(cfg.figDir, 'fig5_convergence.png'));
