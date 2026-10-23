function [alpha, info] = lsc_pcg(y, ctx, dvar, opt)
%LSC_PCG  Matrix-free preconditioned conjugate gradients for (C + D) alpha = y.
%
%   [alpha, info] = LSC_PCG(y, ctx, dvar, opt)
%
%   Solves the same system as LSC_BLOCK_GAUSS_SEIDEL, but with Jacobi
%   preconditioned conjugate gradients instead of a stationary block iteration.
%   Every iteration costs one full covariance matrix-vector product (N^2
%   covariance elements) and the residual decreases monotonically, which makes
%   the convergence trace far more informative than a fixed number of sweeps.
%   It is included as a reference for the discussion in Section 4 of the paper,
%   not as a competitor of the three solvers: the memory is O(N), but the
%   number of matrix-vector products is problem dependent and the iteration is
%   not easily restarted.
%
%   Input
%     y    - N x 1 observations (mGal)
%     ctx  - .N .xyz .r .Gcov
%     dvar - observation noise variance (mGal^2), scalar or N x 1
%     opt  - .maxIter (20) .tol (1e-6, relative residual) .blockSize (4000)
%            .verbose (true) .onIter (callback(alpha, k, relres))
%
%   Output
%     alpha - N x 1 coefficients
%     info  - .iter .time .relres .stop

if nargin < 4, opt = struct(); end
if ~isfield(opt, 'maxIter')  || isempty(opt.maxIter),  opt.maxIter  = 20;   end
if ~isfield(opt, 'tol')      || isempty(opt.tol),      opt.tol      = 1e-6; end
if ~isfield(opt, 'blockSize')|| isempty(opt.blockSize), opt.blockSize = 4000; end
if ~isfield(opt, 'verbose')  || isempty(opt.verbose),  opt.verbose  = true; end
if ~isfield(opt, 'onIter')   || isempty(opt.onIter),   opt.onIter   = [];   end

N = ctx.N;
y = y(:);
if isscalar(dvar), dvar = dvar * ones(N, 1); end
dvar = dvar(:);

t0 = tic;
diagA = cov_diag(ctx) + dvar;
diagA(diagA <= 0) = 1;

x = zeros(N, 1);
r = y;
z = r ./ diagA;
p = z;
rho = r.' * z;
bnorm = norm(y);
relHist = zeros(opt.maxIter, 1);
k = 0;
while k < opt.maxIter
    k   = k + 1;
    Ap  = matvec(p, ctx, dvar, opt.blockSize);
    pAp = p.' * Ap;
    if pAp <= 0, break; end
    a = rho / pAp;
    x = x + a * p;
    r = r - a * Ap;
    rel = norm(r) / bnorm;
    relHist(k) = rel;
    if opt.verbose
        fprintf('   CG %2d  relres = %.3e  t = %.0f s\n', k, rel, toc(t0));
    end
    if ~isempty(opt.onIter), opt.onIter(x, k, rel); end
    if rel < opt.tol, break; end
    z = r ./ diagA;
    rhoNew = r.' * z;
    p = z + (rhoNew / rho) * p;
    rho = rhoNew;
end

alpha = x;
info = struct('iter', k, 'time', toc(t0), 'relres', relHist(1:k), 'stop', 'maxIter');
if k > 0 && relHist(k) < opt.tol, info.stop = 'tol'; end
end

% ----------------------------------------------------------------------
function Ap = matvec(p, ctx, dvar, bs)
N = ctx.N;
Ap = zeros(N, 1);
for i0 = 1:bs:N
    i1 = min(i0 + bs - 1, N);
    ii = (i0:i1).';
    acc = zeros(i1 - i0 + 1, 1);
    for j0 = 1:bs:N
        j1 = min(j0 + bs - 1, N);
        jj = (j0:j1).';
        Cij = tr_cov_pairs(ctx.Gcov, ctx.xyz(ii, :), ctx.r(ii), ...
                           ctx.xyz(jj, :), ctx.r(jj));
        acc = acc + Cij * p(jj);
    end
    Ap(ii) = acc;
end
Ap = Ap + dvar .* p;
end

% ----------------------------------------------------------------------
function d = cov_diag(ctx)
% Diagonal of the covariance matrix, evaluated block by block.
n  = ctx.N;
d  = zeros(n, 1);
bs = 4000;
for i0 = 1:bs:n
    i1 = min(i0 + bs - 1, n);
    ii = (i0:i1).';
    Cii = tr_cov_pairs(ctx.Gcov, ctx.xyz(ii, :), ctx.r(ii), ...
                       ctx.xyz(ii, :), ctx.r(ii));
    d(ii) = diag(Cii);
end
end
