%% MAKE_FIG4_GSVS17  Fig. 4: height-anomaly errors along the GSVS17 profile.
%
%   Reproduces Fig. 4 of the paper: the height-anomaly error of each solver at
%   the 222 GNSS/levelling check points, plotted against the point index, with
%   the ellipsoidal height of the check points on the right-hand axis.
%
%   The datum offset is removed from each error series (it is a property of the
%   reference models, not of the collocation solution) and reported in the
%   command window.
%
%   Solver results are read from results/standard_lsc.mat,
%   results/sequential_lsc.mat, results/block_gs_lsc_2000.mat and
%   results/block_gs_lsc_4000.mat; the ones that are missing are skipped with a
%   message, so the figure can be built incrementally.
%
%   Writes figures/fig4_gsvs17.png and .pdf.
%
%   See also PREDICT_GNSS, RUN_STANDARD_LSC, RUN_SEQUENTIAL_LSC,
%   RUN_BLOCK_GS_LSC.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Fig. 4: height-anomaly errors along GSVS17 ===\n');

D    = load_colorado(cfg);
idx  = (1:numel(D.gnssObs)).';
hgt  = D.gnssBLh(:, 3);            % ellipsoidal height of the check points

% style of each curve, as in the paper: MATLAB default colours.  The curves of
% the three solvers coincide, so decreasing line widths let all of them remain
% visible (the widest, standard LSC, is drawn first).
files = { 'standard_lsc.mat',      'standard',              [0.0000 0.4470 0.7410], 2.6;
          'sequential_lsc.mat',    'sequential',            [0.6350 0.0780 0.1840], 1.8;
          'block_gs_lsc_2000.mat', 'blockGS 2000',          [0.4660 0.6740 0.1880], 1.3;
          'block_gs_lsc_4000.mat', 'blockGS 4000',          [0.9290 0.6940 0.1250], 0.8 };

fig  = figure('Color', 'w', 'Position', [80 80 1000 560]);
hs   = gobjects(0);
labs = {};
yyaxis left
hold on;
for k = 1:size(files, 1)
    fname = fullfile(cfg.resDir, files{k, 1});
    if ~exist(fname, 'file')
        fprintf('   %s not found, skipped\n', files{k, 1});
        continue;
    end
    S = load(fname, 'R');
    R = S.R;

    % rebuild the context used by the run (same observations, same tables)
    ctx = make_context(D.BLh(R.sel, :), cfg, struct('rRef', D.r(R.sel)));
    [~, st] = predict_gnss(ctx.Gcross, R.alpha, ctx, D);

    hs(end + 1)   = plot(idx, st.err * 100, '-', 'LineWidth', files{k, 4}, ...
        'Color', files{k, 3});                                       %#ok<SAGROW>
    labs{end + 1} = files{k, 2};                                     %#ok<SAGROW>
    fprintf('   %-16s STD = %.4f cm   datum = %+.3f m\n', ...
        files{k, 2}, st.std * 100, st.datum);
end
xlim([0 max(idx)]);
ylim([-8 6]);
ylabel('Height anomaly error (cm)');
set(gca, 'FontSize', 12);

yyaxis right
hH = plot(idx, hgt, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
ylim([1800 3400]);
ylabel('Ellipsoidal height (m)');

xlabel('GNSS/levelling point index');
grid on; box on;
ax = gca;
ax.YAxis(1).Color = 'k';           % the paper uses black axes on both sides
ax.YAxis(2).Color = 'k';
ax.XAxis.Color = 'k';

if ~isempty(hs)
    legend([hs, hH], [labs, {'Height'}], 'Orientation', 'horizontal', ...
        'Location', 'northoutside', 'FontSize', 10, 'Box', 'off');
end

print(fig, fullfile(cfg.figDir, 'fig4_gsvs17.png'), '-dpng', '-r300');
print(fig, fullfile(cfg.figDir, 'fig4_gsvs17.pdf'), '-dpdf', '-bestfit');
fprintf('   written %s\n', fullfile(cfg.figDir, 'fig4_gsvs17.png'));
