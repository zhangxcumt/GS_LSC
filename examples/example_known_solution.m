%% EXAMPLE_KNOWN_SOLUTION  Algebraic check of the three collocation solvers.
%
%   This is the quick test of the package.  It needs no input data and runs in
%   a few seconds.
%
%   A synthetic observation vector is built from a known coefficient vector,
%
%       y = (C + D) alpha_true ,
%
%   so that the exact solution of the collocation equations is known.  The
%   three solvers are then applied to y and the recovered coefficients are
%   compared with alpha_true and with each other:
%
%       lsc_standard             exact, Cholesky factorisation of A = C + D
%       lsc_sequential           exact, Moritz/Tscherning recursion
%       lsc_block_gauss_seidel   iterative, no recursion state
%
%   Expected: the two exact solvers reproduce alpha_true to about 1e-12
%   relative, and the block Gauss-Seidel iteration converges to the same
%   vector to about 1e-8.
%
%   Run from the repository root:
%
%       cd examples
%       example_known_solution
%
%   See also GSLSC_CONFIG, MAKE_CONTEXT, LSC_MATVEC.

clearvars -except t0; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Quick test: known-solution check of the three solvers ===\n');

% ----------------------------------------------------------------------
% 1. Synthetic observation geometry (400 points in a 2 deg x 2 deg cell)
% ----------------------------------------------------------------------
rng(7, 'twister');
N   = 400;
lat = 35 + 2 * rand(N, 1);
lon = -106 + 2 * rand(N, 1);
BLh = [lat, lon, zeros(N, 1)];

ctx  = make_context(BLh, cfg);
dvar = cfg.dvar;                       % 16 mGal^2

% ----------------------------------------------------------------------
% 2. Observations from a known coefficient vector
% ----------------------------------------------------------------------
alphaTrue = 0.3 * randn(N, 1);
y = lsc_matvec(alphaTrue, ctx, dvar);
% No mean is removed here: the identity y = (C + D) alpha_true would no longer
% hold, and the solvers themselves do not require a zero-mean right-hand side.
% (Mean removal is a modelling decision taken by the experiment scripts.)

fprintf('   N = %d, D = %g mGal^2, ||y|| = %.4g mGal\n', N, dvar, norm(y));

% ----------------------------------------------------------------------
% 3. Standard LSC
% ----------------------------------------------------------------------
fprintf('\n   [1/3] standard LSC\n');
[aStd, iStd] = lsc_standard(y, ctx, dvar, struct('precision', 'double', ...
    'verbose', false, 'timeMatvec', false));

% ----------------------------------------------------------------------
% 4. Sequential LSC
% ----------------------------------------------------------------------
fprintf('   [2/3] sequential LSC\n');
[aSeq, ~, iSeq] = lsc_sequential(y, ctx, dvar, struct('blockSize', 40, ...
    'single', false, 'verbose', false));

% ----------------------------------------------------------------------
% 5. Block Gauss-Seidel LSC
% ----------------------------------------------------------------------
fprintf('   [3/3] block Gauss-Seidel LSC\n');
[aBgs, ~, iBgs] = lsc_block_gauss_seidel(y, ctx, dvar, struct( ...
    'blockSize', 40, 'stopRule', 'coef', 'tolCoef', 1e-12, 'nSweep', 2000, ...
    'verbose', false));

% ----------------------------------------------------------------------
% 6. Comparison
% ----------------------------------------------------------------------
eStd = norm(aStd - alphaTrue) / norm(alphaTrue);
eSeq = norm(aSeq - alphaTrue) / norm(alphaTrue);
eBgs = norm(aBgs - alphaTrue) / norm(alphaTrue);
dSs  = norm(aSeq - aStd) / norm(aStd);
dBs  = norm(aBgs - aStd) / norm(aStd);

fprintf('\n   %-34s %-12s %-12s\n', 'quantity', 'value', 'tolerance');
fprintf('   %-34s %-12.3e %-12.1e\n', 'rel. error of standard LSC',    eStd, 1e-10);
fprintf('   %-34s %-12.3e %-12.1e\n', 'rel. error of sequential LSC',  eSeq, 1e-10);
fprintf('   %-34s %-12.3e %-12.1e\n', 'rel. error of block GS LSC',    eBgs, 1e-6);
fprintf('   %-34s %-12.3e %-12.1e\n', 'sequential vs standard',        dSs,  1e-11);
fprintf('   %-34s %-12.3e %-12.1e\n', 'block GS vs standard',          dBs,  1e-6);

fprintf('\n   sweeps of the block GS iteration : %d (%s)\n', ...
    iBgs.sweeps, iBgs.stopReason);
fprintf('   exact relative residual          : %.3e\n', iBgs.relResid);
fprintf('   storage  standard : %.2f MB   sequential : %.2f MB   block GS : %.2f MB\n', ...
    iStd.bytes / 2^20, iSeq.pBytes / 2^20, iBgs.blockBytes / 2^20);
fprintf('   time     standard : %.2f s    sequential : %.2f s    block GS : %.2f s\n', ...
    iStd.time, iSeq.time, iBgs.time);

% ----------------------------------------------------------------------
% 7. Verdict
% ----------------------------------------------------------------------
ok = eStd < 1e-10 && eSeq < 1e-10 && dSs < 1e-11 && eBgs < 1e-6 && dBs < 1e-6;
if ok
    fprintf('\n   PASS: the three solvers reproduce the known solution.\n\n');
else
    fprintf('\n   FAIL: see the tolerances above.\n\n');
    error('example_known_solution:failed', 'The solvers do not agree.');
end
