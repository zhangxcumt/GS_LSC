function [alpha, info] = lsc_standard(y, ctx, dvar, opt)
%LSC_STANDARD  Standard (batch) least-squares collocation.
%
%   [alpha, info] = LSC_STANDARD(y, ctx, dvar, opt)
%
%   Solves the collocation normal equations
%
%       (C + D) alpha = y ,        s(P) = sum_j alpha_j C_{sg}(P, Q_j)
%
%   by Cholesky factorisation of the full N x N matrix A = C + D,
%
%       A = R'R ,      alpha = R \ (R' \ y) .
%
%   Cost: O(N^2) memory for A and O(N^3) for the factorisation.  The matrix is
%   stored in single precision by default, which is what allows N = 5 x 10^4 to
%   be factorised on a 32 GB machine (9.3 GB instead of 18.6 GB); the resulting
%   coefficients are accurate to about 1e-6 relative, far below the 1e-2
%   accuracy needed for a height anomaly in cm.
%
%   Input
%     y    - N x 1 observations (residual gravity, mGal), mean removed
%     ctx  - .N .xyz (N x 3) .r (N x 1) .Gcov (TR_COV_GRID of the auto-covariance)
%     dvar - observation noise variance (mGal^2), scalar or N x 1
%     opt  - .precision ('single'|'double', default single)
%            .blockSize (2000)   block size of the matrix assembly
%            .verbose   (true)
%            .timeMatvec (true)  time one matrix-vector product with A before
%                                the factorisation (used for Table 4 of the
%                                paper: the cost of one "cached" sweep)
%            .keepMatrix (false) return A in info.A instead of clearing it
%
%   Output
%     alpha - N x 1 collocation coefficients
%     info  - .precision .buildTime .cholTime .solveTime .time .bytes
%             .nPairEval .tMatvec
%
%   See also LSC_SEQUENTIAL, LSC_BLOCK_GAUSS_SEIDEL, BUILD_COVARIANCE_MATRIX.

if nargin < 4, opt = struct(); end
if ~isfield(opt, 'precision')  || isempty(opt.precision),  opt.precision  = 'single'; end
if ~isfield(opt, 'blockSize')  || isempty(opt.blockSize),  opt.blockSize  = 2000; end
if ~isfield(opt, 'verbose')    || isempty(opt.verbose),    opt.verbose    = true; end
if ~isfield(opt, 'timeMatvec') || isempty(opt.timeMatvec), opt.timeMatvec = true; end
if ~isfield(opt, 'keepMatrix') || isempty(opt.keepMatrix), opt.keepMatrix = false; end

isSingle = strncmpi(opt.precision, 'single', 3);
N = ctx.N;
y = y(:);
if isscalar(dvar), dvar = dvar * ones(N, 1); end

% ---- assemble A = C + D --------------------------------------------------
if opt.verbose, fprintf('   standard LSC: N = %d, assembling A (%s) ...\n', N, opt.precision); end
[A, ainfo] = build_covariance_matrix(ctx, dvar, ...
    struct('blockSize', opt.blockSize, 'single', isSingle, 'verbose', opt.verbose));

if isSingle, y = single(y); end

% ---- cost of one matrix-vector product (Table 4) --------------------------
tMatvec = NaN;
if opt.timeMatvec
    p = single(randn(N, 1));
    acc = 0;                                   %#ok<NASGU>
    t0 = tic;  acc = A * p;  tMatvec = toc(t0); %#ok<NASGU>
    if opt.verbose
        fprintf('   standard LSC: one matvec with A cached = %.3f s\n', tMatvec);
    end
end

% ---- Cholesky factorisation ----------------------------------------------
t0 = tic;
R  = chol(A);
tChol = toc(t0);
clear A
if opt.verbose
    fprintf('   standard LSC: Cholesky %.1f s (%.2f GB)\n', tChol, N * N * 4 / 2^30);
end

% ---- two triangular solves ----------------------------------------------
t0 = tic;
alpha = R \ (R' \ y);
tSolve = toc(t0);
alpha = double(alpha);

info = struct('precision', opt.precision, 'buildTime', ainfo.time, ...
              'cholTime', tChol, 'solveTime', tSolve, ...
              'time', ainfo.time + tChol + tSolve, ...
              'bytes', N * N * 4, 'nPairEval', ainfo.nPairEval, 'tMatvec', tMatvec);
if opt.verbose
    fprintf('   standard LSC: total %.1f s\n', info.time);
end
if opt.keepMatrix, info.R = R; end
end
