function Ap = lsc_matvec(alpha, ctx, dvar, blockSize)
%LSC_MATVEC  Matrix-vector product with the normal-equation matrix A = C + D.
%
%   Ap = LSC_MATVEC(alpha, ctx, dvar, blockSize)
%
%   Evaluates (C + D) * alpha block by block, without ever forming C.  One call
%   costs N^2 covariance evaluations, i.e. exactly one sweep of the block
%   Gauss-Seidel iteration.  It is used
%
%     * to build synthetic observations y = (C + D) alpha_true in the examples,
%     * to measure the exact relative residual ||A alpha - y|| / ||y||.
%
%   Input
%     alpha     - N x 1 vector
%     ctx       - .N .xyz .r .Gcov
%     dvar      - observation noise variance, scalar or N x 1 (mGal^2)
%     blockSize - (optional) block size, default 2000
%
%   Output
%     Ap        - N x 1
%
%   See also LSC_BLOCK_GAUSS_SEIDEL, LSC_PCG.

if nargin < 4 || isempty(blockSize), blockSize = 2000; end

N = ctx.N;
alpha = alpha(:);
if isscalar(dvar), dvar = dvar * ones(N, 1); end
dvar = dvar(:);

Ap = zeros(N, 1, 'like', alpha);
if ~isa(alpha, 'single'), Ap = zeros(N, 1); end

for i0 = 1:blockSize:N
    i1 = min(i0 + blockSize - 1, N);
    ii = (i0:i1).';
    acc = zeros(i1 - i0 + 1, 1, class(alpha));
    for j0 = 1:blockSize:N
        j1 = min(j0 + blockSize - 1, N);
        jj = (j0:j1).';
        Cij = tr_cov_pairs(ctx.Gcov, ctx.xyz(ii, :), ctx.r(ii), ...
                           ctx.xyz(jj, :), ctx.r(jj));
        acc = acc + cast(Cij, class(alpha)) * alpha(jj);
    end
    Ap(ii) = acc;
end
Ap = Ap + cast(dvar, class(alpha)) .* alpha;
end
