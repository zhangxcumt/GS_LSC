function bytes = storage_model(method, N, nb)
%STORAGE_MODEL  Analytic storage requirement of the three collocation solvers.
%
%   bytes = STORAGE_MODEL(method, N, nb)
%
%   Reproduces Fig. 1 of the paper.  The three solvers scale as
%
%       'standard'    4 N^2 bytes   normal-equation matrix A (single precision)
%       'sequential'  8 N^2 bytes   recursion state P (double precision)
%       'blockgs'     8 nb^2 bytes  one diagonal block (double precision)
%
%   The block Gauss-Seidel solver additionally holds O(N) vectors (alpha, y),
%   which are negligible at the scales of interest and are omitted.  For
%   N = 5 x 10^4 the values are 9.31 GB, 18.63 GB and 32 MB (nb = 2000) or
%   128 MB (nb = 4000), as reported in Table 3.
%
%   Input
%     method - 'standard' | 'sequential' | 'blockgs'
%     N      - number of observations (scalar or vector)
%     nb     - block size, used by 'blockgs' only
%
%   Output
%     bytes  - same shape as N
%
%   See also LSC_STANDARD, LSC_SEQUENTIAL, LSC_BLOCK_GAUSS_SEIDEL.

if nargin < 3, nb = 2000; end
switch lower(method)
    case 'standard'
        bytes = 4 * N .* N;
    case 'sequential'
        bytes = 8 * N .* N;
    case {'blockgs', 'block_gauss_seidel'}
        bytes = 8 * nb * nb * ones(size(N));
    otherwise
        error('storage_model:method', 'Unknown method ''%s''.', method);
end
end
