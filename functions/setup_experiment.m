function E = setup_experiment(cfg, nMax)
%SETUP_EXPERIMENT  Common input of every experiment script.
%
%   E = SETUP_EXPERIMENT(cfg) loads the Colorado test field, selects the
%   sub-region of Table 2 of the paper and tabulates the covariance functions.
%
%   E = SETUP_EXPERIMENT(cfg, nMax) keeps only nMax observations of the box
%   (fixed seed, see SELECT_SUBREGION).  nMax = [] keeps all of them.
%
%   The environment variable GSLSC_NMAX overrides nMax.  It is meant for
%   checking the pipeline quickly, e.g.
%
%       set GSLSC_NMAX=20000        (Windows cmd)
%       $env:GSLSC_NMAX=20000       (PowerShell)
%       export GSLSC_NMAX=20000     (bash)
%
%   which runs every experiment on 20 000 observations instead of 50 902.
%
%   Output
%     E.D      structure of LOAD_COLORADO (the whole test field)
%     E.sel    indices of the observations used
%     E.ctx    solver context (covariance tables included)
%     E.y      observations of the sub-region, mean removed (mGal)
%     E.ym     the mean that was removed (mGal)
%     E.dvar   observation noise variance (mGal^2)
%     E.nMax   number of observations requested ([] = all)
%
%   See also GSLSC_CONFIG, LOAD_COLORADO, SELECT_SUBREGION, MAKE_CONTEXT.

if nargin < 2, nMax = []; end
if isempty(nMax), nMax = cfg.nMax; end

env = getenv('GSLSC_NMAX');
if ~isempty(env)
    nMax = str2double(env);
    fprintf('   [GSLSC_NMAX = %d] reduced run, set by the environment variable\n', nMax);
end

fprintf('   loading the Colorado test field ...\n');
D   = load_colorado(cfg);
sel = select_subregion(D, cfg, nMax);

fprintf('   sub-region: %d observations, %.3f-%.3f N / %.3f-%.3f E\n', ...
    numel(sel), min(D.BLh(sel, 1)), max(D.BLh(sel, 1)), ...
    min(D.BLh(sel, 2)), max(D.BLh(sel, 2)));

ctx = make_context(D.BLh(sel, :), cfg, struct('rRef', D.r(sel), 'verbose', true));

yAll = D.y(sel);
ym   = mean(yAll);
y    = yAll - ym;

fprintf('   N = %d, D = %g mGal^2, std(y) = %.4f mGal, mean(y) = %.4f mGal\n', ...
    ctx.N, cfg.dvar, std(yAll), ym);

E = struct('D', D, 'sel', sel, 'ctx', ctx, 'y', y, 'ym', ym, ...
           'dvar', cfg.dvar, 'nMax', nMax);
end
