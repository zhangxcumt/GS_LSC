function [val, T] = signal_predict(Gcross, coef, ctx, xyzP, rP, opt)
%SIGNAL_PREDICT  Predict the signal at arbitrary points from the coefficients.
%
%   T(P)    = sum_j coef_j * C_{T,g}(P, x_j)   disturbing potential
%   zeta(P) = T(P) / gamma(P)                 height anomaly (m)
%
%   The computation is blocked to bound the memory: each block holds
%   opt.pixSize prediction points times opt.chunkSize observations.
%
%   Input
%     Gcross  - cross-covariance interpolant (functional 12 or 13)
%     coef    - N x 1 collocation coefficients
%     ctx     - .N .u .r .xyz
%     xyzP,rP - cartesian coordinates and radii of the prediction points
%     opt     - .pixSize (2000) .chunkSize (4000)
%               .gnorm   normal gravity at the prediction points (mGal); when
%                        given, the height anomaly is returned instead of T
%               .sparse  (true) skip coefficients that are exactly zero
%
%   Output
%     val - height anomaly (m) if opt.gnorm is given, otherwise T
%     T   - disturbing potential

if nargin < 6, opt = struct(); end
if ~isfield(opt,'pixSize')   || isempty(opt.pixSize),   opt.pixSize   = 2000; end
if ~isfield(opt,'chunkSize') || isempty(opt.chunkSize), opt.chunkSize = 4000; end
if ~isfield(opt,'sparse')    || isempty(opt.sparse),    opt.sparse    = true; end

N = ctx.N;
nP = size(xyzP, 1);
T = zeros(nP, 1);
if isfield(ctx,'xyz') && ~isempty(ctx.xyz)
    xyzQ = ctx.xyz;
else
    xyzQ = ctx.u .* ctx.r;          % fallback: rebuild cartesian coordinates
end
rQ = ctx.r;
coef = coef(:);
% Sparse model: only the non-zero coefficients contribute.  For the sparse
% solvers most coefficients are exactly zero, so this reduces the number of
% kernel evaluations of the prediction from N to the number of active ones,
% which matters at large N.
if opt.sparse
    nz = find(coef ~= 0);
    if numel(nz) < N
        coef = coef(nz); xyzQ = xyzQ(nz,:); rQ = rQ(nz); N = numel(nz);
    end
end
i0 = 1;
while i0 <= nP
    i1 = min(i0 + opt.pixSize - 1, nP);
    acc = zeros(i1 - i0 + 1, 1);
    j0 = 1;
    while j0 <= N
        j1 = min(j0 + opt.chunkSize - 1, N);
        C = tr_cov_pairs(Gcross, xyzP(i0:i1,:), rP(i0:i1), xyzQ(j0:j1,:), rQ(j0:j1));
        acc = acc + C * coef(j0:j1);
        j0 = j1 + 1;
    end
    T(i0:i1) = acc;
    i0 = i1 + 1;
end

if isfield(opt,'gnorm') && ~isempty(opt.gnorm)
    val = T ./ opt.gnorm(:);
else
    val = T;
end
end
