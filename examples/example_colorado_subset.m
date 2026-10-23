%% EXAMPLE_COLORADO_SUBSET  The three solvers on a real subsample of Colorado.
%
%   Runs the three solvers on a spatially compact subsample of the Colorado
%   test field (about 6 000 observations, no subsampling beyond the box) and
%   compares the height anomalies predicted at the 222 GNSS/levelling check
%   points.  The purpose is to show that, on real data and with the covariance
%   model of the paper, the three solvers give the same prediction.
%
%   The absolute accuracy of this example is not the accuracy of the paper:
%   it uses about 12 % of the observations of the paper's sub-region, so the
%   residual field is much less constrained.  What is compared is the
%   DIFFERENCE between the solvers, which must vanish.
%
%   Run from the repository root:
%
%       cd examples
%       example_colorado_subset
%
%   See also SELECT_SUBREGION, LSC_STANDARD, LSC_SEQUENTIAL,
%   LSC_BLOCK_GAUSS_SEIDEL, PREDICT_GNSS.

clearvars -except t0; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));
cfg = gslsc_config();
fprintf('\n=== Colorado subsample: the three solvers on real data ===\n');

% ----------------------------------------------------------------------
% 1. Data and a compact subsample
% ----------------------------------------------------------------------
fprintf('   loading the Colorado test field ...\n');
D = load_colorado(cfg);
fprintf('   %d observations in the test field, %d GNSS/levelling check points\n', ...
    D.n, numel(D.gnssObs));

cfgSub    = cfg;
cfgSub.nLat = [37.60 38.80];
cfgSub.nLon = [-106.80 -105.60];
sel = select_subregion(D, cfgSub);
fprintf('   %d observations inside %.2f-%.2f N / %.2f-%.2f E\n', ...
    numel(sel), cfgSub.nLat, cfgSub.nLon);

ctx  = make_context(D.BLh(sel, :), cfg, struct('rRef', D.r));
dvar = cfg.dvar;
y    = D.y(sel);
y    = y - mean(y);
fprintf('   N = %d, std(y) = %.3f mGal\n', ctx.N, std(y));

nb = 400;                                   % block size of the two block solvers

% The stopping rule of the paper: the iteration ends when the largest change of
% the predicted height anomaly at the check points falls below 0.01 cm.
xyzG = BLh_XYZR(D.gnssBLh);
thetaFun = @(a) 100 * signal_predict(ctx.Gcross, a, ctx, ...
    xyzG(:, 1:3), xyzG(:, 4), struct('gnorm', D.gnssGnorm));

% ----------------------------------------------------------------------
% 2. The three solutions
% ----------------------------------------------------------------------
fprintf('\n   [1/3] standard LSC\n');
[aStd, iStd] = lsc_standard(y, ctx, dvar, struct('precision', 'double', ...
    'blockSize', nb, 'verbose', false, 'timeMatvec', false));

fprintf('   [2/3] sequential LSC\n');
[aSeq, ~, iSeq] = lsc_sequential(y, ctx, dvar, struct('blockSize', nb, ...
    'single', false, 'verbose', false));

fprintf('   [3/3] block Gauss-Seidel LSC (stops on the predicted quantity)\n');
[aBgs, ~, iBgs] = lsc_block_gauss_seidel(y, ctx, dvar, struct( ...
    'blockSize', nb, 'stopRule', 'pred', 'tolPred', 0.01, 'thetaFun', thetaFun, ...
    'nSweep', 500, 'verbose', false));

% ----------------------------------------------------------------------
% 3. Prediction at the GNSS/levelling check points
% ----------------------------------------------------------------------
[~, SL] = predict_gnss(ctx.Gcross, aStd, ctx, D);
[~, SS] = predict_gnss(ctx.Gcross, aSeq, ctx, D);
[~, SB] = predict_gnss(ctx.Gcross, aBgs, ctx, D);

fprintf('\n   %-30s %-12s %-12s\n', 'solver', 'STD (cm)', 'time (s)');
fprintf('   %-30s %-12.4f %-12.2f\n', 'standard LSC',            SL.std * 100, iStd.time);
fprintf('   %-30s %-12.4f %-12.2f\n', 'sequential LSC',          SS.std * 100, iSeq.time);
fprintf('   %-30s %-12.4f %-12.2f\n', 'block Gauss-Seidel LSC',  SB.std * 100, iBgs.time);

dSeq = max(abs(SS.err - SL.err)) * 100;      % cm
dBgs = max(abs(SB.err - SL.err)) * 100;      % cm
dStd = abs(SB.std - SL.std) * 100;           % cm
fprintf('\n   max |dzeta|, sequential vs standard      : %.3e cm\n', dSeq);
fprintf('   max |dzeta|, block GS  vs standard      : %.3e cm\n', dBgs);
fprintf('   |STD(block GS) - STD(standard)|         : %.3e cm\n', dStd);
fprintf('   sweeps of the block GS iteration        : %d (%s)\n', ...
    iBgs.sweeps, iBgs.stopReason);
fprintf('   exact relative residual of block GS     : %.3e\n', iBgs.relResid);

% The block Gauss-Seidel solution is stopped when the height anomaly at the
% check points changes by less than 0.01 cm between two sweeps.  The remaining
% difference to the exact solution is the last change accumulated over the
% final sweeps, so it is of the order of the threshold and not of it squared;
% this is the behaviour discussed in Section 4.1 of the paper.  The statistic
% that is compared in Table 3 (the standard deviation at the check points) is
% insensitive to it.

% ----------------------------------------------------------------------
% 4. Verdict
% ----------------------------------------------------------------------
ok = dSeq < 1e-6 && dStd < 0.01 && dBgs < 0.05;   % cm
if ok
    fprintf('\n   PASS: the three solvers predict the same height anomaly.\n\n');
else
    fprintf('\n   FAIL: the solvers disagree by more than the tolerance.\n\n');
    error('example_colorado_subset:failed', 'The solvers do not agree.');
end
