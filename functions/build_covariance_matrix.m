function [A, info] = build_covariance_matrix(ctx, dvar, opt)
%BUILD_COVARIANCE_MATRIX  Assemble the full N x N normal-equation matrix A = C + D.
%
%   [A, info] = BUILD_COVARIANCE_MATRIX(ctx, dvar, opt)
%
%   Standard collocation needs the explicit normal-equation matrix.  Assembling
%   it naively creates a double-precision copy of the whole matrix before it is
%   cast, which doubles the peak memory (20 GB instead of 10 GB at N = 5 x 10^4).
%   This function therefore writes the matrix block by block into a
%   pre-allocated array; the temporary is one block row (n_b x N), so the peak
%   memory is the matrix itself.
%
%   Input
%     ctx   - .N .xyz .r .Gcov     (see LSC_STANDARD)
%     dvar  - observation noise variance, scalar or N x 1 (mGal^2)
%     opt   - .blockSize (2000) .single (true) .verbose (true)
%
%   Output
%     A     - N x N symmetric positive definite matrix, single or double
%     info  - .time .bytes .nPairEval .class .blockSize
%
%   See also LSC_STANDARD, TR_COV_PAIRS.

if nargin < 3, opt = struct(); end
if ~isfield(opt, 'blockSize') || isempty(opt.blockSize), opt.blockSize = 2000; end
if ~isfield(opt, 'single')    || isempty(opt.single),    opt.single    = true;  end
if ~isfield(opt, 'verbose')   || isempty(opt.verbose),   opt.verbose   = true;  end

N = ctx.N;
if isscalar(dvar), dvar = dvar * ones(N, 1); end
dvar = dvar(:);
cls  = 'double';
if opt.single, cls = 'single'; end

A = zeros(N, N, cls);
r   = ctx.r(:);
xyz = ctx.xyz;

starts = 1:opt.blockSize:N;
ends   = min(starts + opt.blockSize - 1, N);
nb     = numel(starts);
tStart = tic;
nPair  = 0;

% Only the upper triangle is evaluated; the lower one is a copy of it.
for i = 1:nb
    ik = (starts(i):ends(i)).';
    for j = i:nb
        ij  = (starts(j):ends(j)).';
        Cij = tr_cov_pairs(ctx.Gcov, xyz(ik, :), r(ik), xyz(ij, :), r(ij));
        nPair = nPair + numel(ik) * numel(ij);
        A(ik, ij) = cast(Cij, cls);
        if j > i
            A(ij, ik) = cast(Cij.', cls);
        end
    end
    if opt.verbose && (mod(i, max(1, ceil(nb / 10))) == 0 || i == nb)
        fprintf('   build A  %3d/%3d  t = %6.1f s\n', i, nb, toc(tStart));
    end
end

dg = (1:N+1:N^2).';
A(dg) = A(dg) + cast(dvar, cls);

info = struct('time', toc(tStart), 'bytes', N * N * (4 + 4 * ~opt.single), ...
              'nPairEval', nPair, 'class', cls, 'blockSize', opt.blockSize);
end
