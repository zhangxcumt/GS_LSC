%% RUN_QUICK_TEST  Run the quick test of the package.
%
%   This is the test referred to by the "Computer Code Availability" section of
%   the manuscript.  It runs, in this order,
%
%       1. example_known_solution   algebraic check against a known solution,
%                                   no input data needed (a few seconds)
%       2. example_colorado_subset  the three solvers on ~2 000 real Colorado
%                                   observations (about one minute)
%
%   and stops with an error if any of the two fails.
%
%   Run from the repository root:
%
%       cd examples
%       run_quick_test
%
%   or from the MATLAB command window:
%
%       run('examples/run_quick_test.m')
%
%   Requirements: MATLAB R2018b or newer, no toolbox.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

t0 = tic;
fprintf('\n################ QUICK TEST ################\n');

run(fullfile(fileparts(mfilename('fullpath')), 'example_known_solution.m'));
run(fullfile(fileparts(mfilename('fullpath')), 'example_colorado_subset.m'));

fprintf('################ QUICK TEST PASSED ################\n');
fprintf('total elapsed time: %.1f s\n\n', toc(t0));
