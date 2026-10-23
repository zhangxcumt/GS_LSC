function [alpha, E, info] = lsc_block_gauss_seidel(y, ctx, dvar, opt)
%LSC_BLOCK_GAUSS_SEIDEL  Block Gauss-Seidel solution of the collocation equations.
%
%   [alpha, E, info] = LSC_BLOCK_GAUSS_SEIDEL(y, ctx, dvar, opt)
%
%   Solves (C + D) alpha = y directly, keeping no recursion state:
%
%       alpha_k^(m+1) = A_kk^{-1} ( y_k - sum_{j<k} A_kj alpha_j^(m+1)
%                                        - sum_{j>k} A_kj alpha_j^(m) ) ,
%       A_kj = C_kj + delta_kj D_k .
%
%   Only the n_b x n_b diagonal blocks are factorised, so the storage is
%   O(n_b^2 + N) and independent of N (32 MB for n_b = 2000, 128 MB for
%   n_b = 4000), against O(N^2) for both standard and sequential LSC.  The
%   price is that every sweep evaluates all N^2 covariance elements again: the
%   method trades time for storage, it is not a "faster sequential" method.
%
%   STOPPING.  The coefficients converge much more slowly than the predicted
%   quantity, so stopping on the residual or on the change of alpha overestimates
%   the number of sweeps by more than an order of magnitude.  The default rule
%   therefore stops on the predicted quantity: the iteration ends when the
%   largest change of the height anomaly at the check points between two
%   successive sweeps is below opt.tolPred (0.01 cm in the paper).  Supply the
%   prediction through opt.thetaFun, a handle that maps alpha to the vector of
%   predicted quantities:
%
%       thetaFun = @(a) predict_zeta(a) - zetaRef ;
%
%   Input
%     y    - N x 1 observations (mGal), mean removed
%     ctx  - .N .u .xyz .r .Gcov
%     dvar - observation noise variance (mGal^2), scalar or N x 1
%     opt  - .blockSize (2000)  n_b
%            .order      ('lat')  'lat' | 'tile' | 'asGiven'
%            .stopRule   ('pred') 'pred' | 'coef' | 'sweeps'
%            .tolPred    (0.01)   threshold on max|d theta|, same unit as thetaFun
%            .thetaFun   ([])     handle alpha -> vector of predicted quantities
%            .tolCoef    (1e-6)   relative change of alpha, rule 'coef'
%            .nSweep     (100)    maximum number of sweeps
%            .omega      (1)      relaxation factor (1 = Gauss-Seidel)
%            .single     (false)  single-precision covariance blocks
%            .alpha0     ([])     warm start
%            .onSweep    ([])     callback(alpha, m, dAlpha, theta, dTheta)
%            .residualEvery (0)   evaluate the exact residual every k sweeps
%                                 (costs one matrix-vector product each time;
%                                 0 keeps only the final residual)
%            .verbose    (true)
%
%   Output
%     alpha - N x 1 coefficients, restored to the input observation order
%     E     - N x 1 self-fitting residual, E = model - y = -D.*alpha
%     info  - .time .nBlock .blockSize .nPairEval .sweeps .stopReason
%             .dAlphaHist .residHist .dThetaHist .sweepTimes .meanSweepTime
%             .relResid (exact, at the last sweep) .theta .blockBytes
%
%   See also LSC_STANDARD, LSC_SEQUENTIAL, LSC_PCG.

if nargin < 4, opt = struct(); end
if ~isfield(opt, 'blockSize') || isempty(opt.blockSize), opt.blockSize = 2000; end
if ~isfield(opt, 'order')     || isempty(opt.order),     opt.order     = 'lat';  end
if ~isfield(opt, 'stopRule')  || isempty(opt.stopRule),  opt.stopRule  = 'pred'; end
if ~isfield(opt, 'tolPred')   || isempty(opt.tolPred),   opt.tolPred   = 0.01;   end
if ~isfield(opt, 'thetaFun')  || isempty(opt.thetaFun),  opt.thetaFun  = [];     end
if ~isfield(opt, 'tolCoef')   || isempty(opt.tolCoef),   opt.tolCoef   = 1e-6;   end
if ~isfield(opt, 'nSweep')    || isempty(opt.nSweep),    opt.nSweep    = 100;    end
if ~isfield(opt, 'verbose')   || isempty(opt.verbose),   opt.verbose   = true;   end
if ~isfield(opt, 'single')    || isempty(opt.single),    opt.single    = false;  end
if ~isfield(opt, 'omega')     || isempty(opt.omega),     opt.omega     = 1;      end
if ~isfield(opt, 'onSweep')   || isempty(opt.onSweep),   opt.onSweep   = [];     end
if ~isfield(opt, 'alpha0')    || isempty(opt.alpha0),    opt.alpha0    = [];     end
if ~isfield(opt, 'residualEvery') || isempty(opt.residualEvery), opt.residualEvery = 0; end
cl = @(x) x;
if opt.single, cl = @(x) single(x); end

if strcmpi(opt.stopRule, 'pred') && isempty(opt.thetaFun)
    error('lsc_block_gauss_seidel:thetaFun', ...
        ['stopRule = ''pred'' needs opt.thetaFun, a handle that maps the ' ...
         'coefficient vector to the predicted quantity.\n  Use ' ...
         'stopRule = ''coef'' or ''sweeps'' to stop without prediction.']);
end

N = ctx.N;
y = y(:);
if isscalar(dvar), dvar = dvar * ones(N, 1); end
dvar = dvar(:);

% ---- ordering: 'lat' gives latitude strips, 'tile' gives compact tiles ----
switch lower(opt.order)
    case 'lat'
        lonDeg = atan2d(ctx.u(:, 2), ctx.u(:, 1));
        latDeg = asind(ctx.u(:, 3));
        [~, ord] = sortrows([latDeg, lonDeg]);
    case 'tile'
        lonDeg = atan2d(ctx.u(:, 2), ctx.u(:, 1));
        latDeg = asind(ctx.u(:, 3));
        nSide = max(1, ceil(sqrt(N / opt.blockSize)));
        ix = floor((lonDeg - min(lonDeg)) ./ (max(lonDeg) - min(lonDeg) + eps) * nSide);
        iy = floor((latDeg - min(latDeg)) ./ (max(latDeg) - min(latDeg) + eps) * nSide);
        ix = min(ix, nSide - 1);
        iy = min(iy, nSide - 1);
        [~, ord] = sortrows([iy, ix]);
    otherwise
        ord = (1:N).';
end
back(ord) = 1:N;

if isfield(ctx, 'xyz') && ~isempty(ctx.xyz)
    xyz = ctx.xyz(ord, :);
else
    xyz = ctx.u(ord, :) .* ctx.r(ord);
end
r  = ctx.r(ord);  r  = r(:);
ys = cl(y(ord));     ys = ys(:);
ds = cl(dvar(ord));  ds = ds(:);

starts = 1:opt.blockSize:N;
ends   = min(starts + opt.blockSize - 1, N);
nb     = numel(starts);

alpha = cl(zeros(N, 1));
if ~isempty(opt.alpha0)
    alpha = opt.alpha0(:);
    if numel(alpha) ~= N
        error('lsc_block_gauss_seidel:alpha0', ...
            'alpha0 has %d elements but N = %d.', numel(alpha), N);
    end
    alpha = cl(alpha(ord));
end

tStart = tic;
nPair  = 0;
stopReason = 'maxSweep';
dAlphaHist = nan(opt.nSweep, 1);
residHist  = nan(opt.nSweep, 1);
dThetaHist = nan(opt.nSweep, 1);
sweepTimes = nan(opt.nSweep, 1);
theta = [];
nSweepDone = 0;
tPrev = 0;

for m = 1:opt.nSweep
    aOld = alpha;
    thetaOld = theta;

    for k = 1:nb
        ik = (starts(k):ends(k)).';
        nk = numel(ik);
        s  = cl(zeros(nk, 1));
        for j = 1:nb
            if j == k, continue; end
            ij  = (starts(j):ends(j)).';
            Ckj = cl(tr_cov_pairs(ctx.Gcov, xyz(ik, :), r(ik), xyz(ij, :), r(ij)));
            nPair = nPair + nk * numel(ij);
            s = s + Ckj * alpha(ij);
        end
        Ckk = cl(tr_cov_pairs(ctx.Gcov, xyz(ik, :), r(ik), xyz(ik, :), r(ik)));
        dg  = (1:(nk+1):numel(Ckk)).';
        Ckk(dg) = Ckk(dg) + ds(ik);

        rhs  = ys(ik) - s;
        prev = alpha(ik);
        alpha(ik) = prev + opt.omega * (Ckk \ rhs - prev);
    end

    nSweepDone = m;
    sweepTimes(m) = toc(tStart) - tPrev;
    tPrev = toc(tStart);
    dA = max(abs(alpha - aOld)) / max(max(abs(alpha)), eps);
    dAlphaHist(m) = dA;

    % The residual is only meaningful when it is evaluated for the whole sweep:
    % inside the sweep every block is solved exactly, so a per-block residual is
    % zero by construction.  opt.residualEvery > 0 requests the exact residual
    % every so many sweeps, at the cost of one matrix-vector product each time.
    residHist(m) = NaN;
    if opt.residualEvery > 0 && (mod(m, opt.residualEvery) == 0 || m == opt.nSweep)
        residHist(m) = relresidual(alpha(back), y, ctx, dvar, opt.blockSize);
    end

    % ---- predicted quantity (height anomaly at the check points) ----------
    dTheta = NaN;
    if ~isempty(opt.thetaFun)
        theta  = opt.thetaFun(alpha(back));
        if ~isempty(thetaOld)
            dTheta = max(abs(theta - thetaOld));
        end
        dThetaHist(m) = dTheta;
    end

    if opt.verbose
        fprintf('   BGS sweep %3d  dAlpha = %.3e', m, dA);
        if ~isnan(residHist(m)), fprintf('  resid = %.3e', residHist(m)); end
        if ~isnan(dTheta), fprintf('  max|dtheta| = %.4f', dTheta); end
        fprintf('  t = %6.1f s\n', toc(tStart));
    end
    if ~isempty(opt.onSweep)
        opt.onSweep(alpha(back), m, dA, theta, dTheta);
    end

    % ---- stopping test ----------------------------------------------------
    switch lower(opt.stopRule)
        case 'pred'
            if m > 1 && ~isnan(dTheta) && dTheta < opt.tolPred
                stopReason = 'predicted'; break;
            end
        case 'coef'
            if dA < opt.tolCoef
                stopReason = 'coefficients'; break;
            end
        otherwise
            % fixed number of sweeps
    end
end

alpha = alpha(back);
E = -dvar .* double(alpha);

% ---- exact relative residual (one extra matrix-vector product) -----------
t0 = tic;
res = relresidual(alpha, y, ctx, dvar, opt.blockSize);
residualTime = toc(t0);
if opt.verbose
    fprintf('   BGS: %d sweeps, %s, ||A alpha - y||/||y|| = %.3e, %.1f s\n', ...
        nSweepDone, stopReason, res, toc(tStart));
end

info = struct('time', toc(tStart), 'nBlock', nb, 'blockSize', opt.blockSize, ...
              'nPairEval', nPair, 'sweeps', nSweepDone, 'stopReason', stopReason, ...
              'dAlphaHist', dAlphaHist(1:nSweepDone), ...
              'residHist',  residHist(1:nSweepDone), ...
              'dThetaHist', dThetaHist(1:nSweepDone), ...
              'sweepTimes', sweepTimes(1:nSweepDone), ...
              'meanSweepTime', mean(sweepTimes(1:nSweepDone)), ...
              'residualTime', residualTime, ...
              'relResid', res, 'theta', theta, ...
              'blockBytes', opt.blockSize^2 * 8);
end

% ----------------------------------------------------------------------
function rel = relresidual(alpha, y, ctx, dvar, bs)
% Exact relative residual ||(C + D) alpha - y|| / ||y||.
Ap  = lsc_matvec(alpha, ctx, dvar, bs);
rel = norm(Ap - y(:)) / norm(y);
end
