%% RUN_ALL  Reproduce every table and figure of the paper.
%
%   Runs, in order,
%
%       run_standard_lsc         Table 3, row 1     (~2 min, 20 GB peak)
%       run_sequential_lsc       Table 3, row 2     (~25 min, 20 GB)
%       run_block_gs_lsc         Table 3, rows 3-4  (~35 min, 128 MB)
%       make_table3              Table 3 on screen and in results/table3.csv
%       run_table4_sweep_cost    Table 4
%       make_fig1_storage        Fig. 1
%       make_fig4_gsvs17         Fig. 4
%       make_fig5_convergence    Fig. 5              (~1.5 h)
%
%   The times are those of the 31.4 GB machine of the paper at N = 50 902.
%   Nothing needs a toolbox.  For a quick check of the whole chain, set
%
%       GSLSC_NMAX = 20000
%
%   in the environment: every experiment then runs on 20 000 observations and
%   the whole sequence takes about five minutes.
%
%   See also RUN_QUICK_TEST (examples/).

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
tAll = tic;

steps = {'run_standard_lsc', 'run_sequential_lsc', 'run_block_gs_lsc', ...
         'make_table3', 'run_table4_sweep_cost', 'make_fig1_storage', ...
         'make_fig4_gsvs17', 'make_fig5_convergence'};

for k = 1:numel(steps)
    t0 = tic;
    fprintf('\n############ %s ############\n', steps{k});
    run(fullfile(fileparts(mfilename('fullpath')), [steps{k} '.m']));
    fprintf('############ %s done in %.1f s ############\n', steps{k}, toc(t0));
end

fprintf('\nAll tables and figures produced in %.1f min.\n', toc(tAll) / 60);
fprintf('results/: %s\n', cfg.resDir);
fprintf('figures/: %s\n', cfg.figDir);
