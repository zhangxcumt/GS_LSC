function plot_fig5(store, nbList, cfg)
%PLOT_FIG5  Draw Fig. 5 of the paper from the convergence histories.
%
%   PLOT_FIG5(store, nbList, cfg) plots, for every block size in nbList, the
%   largest change of the predicted height anomaly between two successive
%   sweeps, together with the 0.01 cm stopping threshold.
%
%   store{k} must have the field .dTheta, as produced by MAKE_FIG5_CONVERGENCE
%   or read from the .info.dThetaHist field of the block Gauss-Seidel results.
%   The figure is written to figures/fig5_convergence.png and .pdf.
%
%   See also MAKE_FIG5_CONVERGENCE.

% colours and styles of the paper: n_b = 2000 dark red solid, n_b = 4000 blue
% dashed, tolerance grey dotted
cols  = [0.6350 0.0780 0.1840; 0.0000 0.4470 0.7410; 0.4660 0.6740 0.1880];
style = {'-', '--', '-'};

fig = figure('Color', 'w', 'Position', [100 100 780 520]);
hs  = gobjects(0);
labs = {};
hold on;
for k = 1:numel(nbList)
    yv = store{k}.dTheta;
    if isempty(yv), continue; end
    yv = yv(:);
    yv(yv <= 0) = NaN;                 % log axis: drop the round-off floor
    hs(end + 1) = semilogy(1:numel(yv), yv, style{mod(k - 1, numel(style)) + 1}, ...
        'LineWidth', 2.0, 'Color', cols(mod(k - 1, size(cols, 1)) + 1, :)); %#ok<SAGROW>
    labs{end + 1} = sprintf('n_b = %d', nbList(k));                  %#ok<SAGROW>
end
ht = yline(0.01, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5);

set(gca, 'YScale', 'log', 'FontSize', 12, 'XMinorGrid', 'off', 'YMinorGrid', 'on', ...
    'MinorGridLineStyle', ':', 'GridLineStyle', '-', 'GridAlpha', 0.15);
grid on; box on;

xmax = 0;
for k = 1:numel(store)
    if ~isempty(store{k}.dTheta), xmax = max(xmax, numel(store{k}.dTheta)); end
end
xlim([0, max(40, xmax)]);
ylim([8e-3, 1e1]);
xlabel('Sweep index');
ylabel('Change of predicted height anomaly (cm)');
if ~isempty(hs)
    legend([hs, ht], [labs, {'Tolerance 0.01 cm'}], 'Location', 'northeast', ...
        'FontSize', 10);
end

print(fig, fullfile(cfg.figDir, 'fig5_convergence.png'), '-dpng', '-r300');
print(fig, fullfile(cfg.figDir, 'fig5_convergence.pdf'), '-dpdf', '-bestfit');
end
