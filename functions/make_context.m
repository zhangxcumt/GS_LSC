function ctx = make_context(BLh, cfg, opt)
%MAKE_CONTEXT  Covariance tables and solver context for a set of points.
%
%   ctx = MAKE_CONTEXT(BLh, cfg) builds everything the three solvers need:
%   the tabulated auto-covariance (gravity - gravity) and T - gravity
%   cross-covariance, and the geometry of the observation points.
%
%   Input
%     BLh - n x 3 [latitude longitude height] (deg, deg, m)
%     cfg - structure from GSLSC_CONFIG
%     opt - .rRef     radii used to set the adaptive range of the s grid
%                      (default: the radii of BLh)
%           .psiMax   maximum spherical distance of the table (deg)
%           .verbose  (false)
%
%   Output
%     ctx.N       number of points
%     ctx.u       n x 3 unit vectors
%     ctx.xyz     n x 3 geocentric cartesian coordinates
%     ctx.r       n x 1 geocentric radii
%     ctx.RB      Bjerhammar radius
%     ctx.Gcov    TR_COV_GRID of the auto-covariance (cfg.fun)
%     ctx.Gcross  TR_COV_GRID of the T - gravity cross-covariance (cfg.crossFun)
%
%   See also LOAD_COLORADO, SELECT_SUBREGION.

if nargin < 3, opt = struct(); end
if ~isfield(opt, 'rRef')    || isempty(opt.rRef),    opt.rRef = []; end
if ~isfield(opt, 'psiMax')  || isempty(opt.psiMax),  opt.psiMax = cfg.cov.psiMax; end
if ~isfield(opt, 'verbose') || isempty(opt.verbose), opt.verbose = false; end

xyzr = BLh_XYZR(BLh);
n    = size(BLh, 1);
rRef = opt.rRef;
if isempty(rRef), rRef = xyzr(:, 4); end

P = cfg.cov;
gopt = struct('rRef', rRef, 'psiMax', opt.psiMax);

if opt.verbose, fprintf('   tabulating covariances ...\n'); end
Gcov   = tr_cov_grid(P, cfg.fun,      gopt);
Gcross = tr_cov_grid(P, cfg.crossFun, gopt);
if opt.verbose
    fprintf('   covariance table: resT = %d, resS = %d, psiMax = %g deg\n', ...
        Gcov.info.resT, Gcov.info.resS, Gcov.info.psiMax);
end

ctx = struct();
ctx.N      = n;
ctx.u      = xyzr(:, 1:3) ./ xyzr(:, 4);
ctx.xyz    = xyzr(:, 1:3);
ctx.r      = xyzr(:, 4);
ctx.RB     = P.RB;
ctx.Gcov   = Gcov;
ctx.Gcross = Gcross;
end
