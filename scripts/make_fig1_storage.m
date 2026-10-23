%% MAKE_FIG1_STORAGE  Fig. 1: storage of the three solvers against N.
%
%   Reproduces Fig. 1 of the paper: the storage of the three algorithms as a
%   function of the number of observations (see STORAGE_MODEL),
%
%       standard LSC           4 N^2 bytes   (single-precision matrix A)
%       sequential LSC         8 N^2 bytes   (double-precision state P)
%       block Gauss-Seidel     8 nb^2 bytes  (one diagonal block, nb = 4000)
%
%   The grey dotted line is the 31.4 GB memory of the computer used in the
%   paper.
%
%   Writes figures/fig1_storage.png and .pdf.
%
%   See also STORAGE_MODEL, MAKE_TABLE3.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Fig. 1: storage against N ===\n');

N     = logspace(log10(2.4e3), log10(3.1e5), 500);
GB    = 2^30;
sStd  = storage_model('standard',   N, [])    / GB;
sSeq  = storage_model('sequential', N, [])    / GB;
sBGS  = storage_model('blockgs',    N, 4000)  / GB;

fig = figure('Color', 'w', 'Position', [100 100 760 500]);
loglog(N, sStd, '-',  'LineWidth', 2.0, 'Color', [0.6350 0.0780 0.1840]); hold on;
loglog(N, sSeq, '--', 'LineWidth', 2.0, 'Color', [0.0000 0.4470 0.7410]);
loglog(N, sBGS, '-',  'LineWidth', 2.0, 'Color', [0.4660 0.6740 0.1880]);
yline(31.4, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5);

set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
    'MinorGridLineStyle', ':', 'GridLineStyle', '-', 'GridAlpha', 0.15);
grid on; box on;
xlim([2.4e3 3.1e5]);
ylim([1e-3 1e3]);
xlabel('Number of observations N');
ylabel('Memory (GB)');
legend({'Standard LSC  (A, single)', 'Sequential LSC  (P, double)', ...
        'Block Gauss-Seidel  (n_b = 4000)'}, 'Location', 'northwest', 'FontSize', 10);

print(fig, fullfile(cfg.figDir, 'fig1_storage.png'), '-dpng', '-r300');
print(fig, fullfile(cfg.figDir, 'fig1_storage.pdf'), '-dpdf', '-bestfit');
fprintf('   written %s\n', fullfile(cfg.figDir, 'fig1_storage.png'));
