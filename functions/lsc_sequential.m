function [alpha, E, info] = lsc_sequential(y, ctx, dvar, opt)
%LSC_SEQUENTIAL  Sequential (stepwise) least-squares collocation.
%
%   [alpha, E, info] = LSC_SEQUENTIAL(y, ctx, dvar, opt)
%
%   Moritz's sequential collocation (Moritz 1973, OSU Report 203; Tscherning
%   1974, OSU Report 212).  The observations are sorted into spatially
%   contiguous groups of n_b and processed recursively, so that the only matrix
%   inverts is n_b x n_b:
%
%       v_i    = y_i - A_{i,<i} alpha_{<i}                  innovation
%       S_i    = A_ii - A_{i,<i} P_{<i} A_{<i,i}            n_b x n_b
%       K_i    = P_{<i} A_{<i,i} S_i^{-1}                   gain
%       alpha_{<i} <- alpha_{<i} - K_i v_i
%       alpha_i     = S_i^{-1} v_i
%       P_{<=i} = [ P_{<i} + K_i S_i K_i'   -K_i   ]
%                 [ -K_i'                   S_i^{-1} ]
%
%   with A = C + D and P_{<i} = (A_{<i,<i})^{-1} the inverse covariance of the
%   observations processed so far.  The sign of the gain update is a minus
%   sign, because the state is the coefficient vector A^{-1}y and not a random
%   signal (opposite to the Kalman filter).
%
%   STORAGE.  The matrix to be inverted is only n_b x n_b, but the recursion
%   state P is N x N: sequential collocation does not reduce the total storage
%   requirement.  This is stated by Tscherning (1974): "the total storage
%   requirements are not diminished ... stepwise solution of the normal
%   equations is of no real advantage", and it is quantified in the paper.  The
%   function refuses to run when P would exceed opt.maxMemGB.
%
%   The solution is exact after one pass: no iteration and no convergence
%   criterion are involved.
%
%   MEMORY NOTE.  The gain update of the upper-left block is
%
%       P + K S K'  =  P + (P A) S^{-1} (P A)' ,
%
%   written here as (PA*Si)*PA' and accumulated row block by row block.  The
%   naive form K*S*K' materialises a prev x prev temporary (16.9 GB in double
%   precision at N = 5 x 10^4), which together with P exceeds the memory of the
%   machine; the blocked form needs only bw x prev.
%
%   Input
%     y    - N x 1 observations (mGal), mean removed
%     ctx  - .N .xyz .r .Gcov
%     dvar - observation noise variance (mGal^2), scalar or N x 1
%     opt  - .blockSize (2000)  n_b
%            .order     ('lat') spatial ordering, 'lat' or 'asGiven'
%            .single    (false) store and update P in single precision
%            .maxMemGB  (26)   refuse to run if P needs more than this
%            .verbose   (true)
%
%   Output
%     alpha - N x 1 coefficients, restored to the input observation order
%     E     - N x 1 self-fitting residual, E = model - y = -D.*alpha
%     info  - .time .nBlock .blockSize .nPairEval .pBytes .single .exact
%
%   See also LSC_STANDARD, LSC_BLOCK_GAUSS_SEIDEL.

if nargin < 4, opt = struct(); end
if ~isfield(opt, 'blockSize') || isempty(opt.blockSize), opt.blockSize = 2000; end
if ~isfield(opt, 'order')     || isempty(opt.order),     opt.order     = 'lat';  end
if ~isfield(opt, 'verbose')   || isempty(opt.verbose),   opt.verbose   = true;  end
if ~isfield(opt, 'maxMemGB')  || isempty(opt.maxMemGB),  opt.maxMemGB  = 26;    end
if ~isfield(opt, 'single')    || isempty(opt.single),    opt.single    = false; end
cl = @(x) x;
if opt.single, cl = @(x) single(x); end

N = ctx.N;
y = y(:);
if isscalar(dvar), dvar = dvar * ones(N, 1); end
dvar = dvar(:);

% ---- storage of the recursion state --------------------------------------
pBytes = N * N * 8;
if opt.single, pBytes = N * N * 4; end
if pBytes / 1024^3 > opt.maxMemGB
    error('lsc_sequential:mem', ...
        ['The recursion state P needs %.1f GB (N = %d), more than the %.1f GB ' ...
         'allowed.\n  Sequential LSC has an O(N^2) state that cannot be ' ...
         'compressed (Tscherning 1974).\n  Run a smaller subset, or use ' ...
         'LSC_BLOCK_GAUSS_SEIDEL, whose storage is O(nb^2 + N).'], ...
        pBytes / 1024^3, N, opt.maxMemGB);
end
if opt.verbose
    fprintf('   sequential LSC: N = %d, recursion state P = %.2f GB\n', N, pBytes / 1024^3);
end

% ---- spatial ordering (identical to LSC_BLOCK_GAUSS_SEIDEL) --------------
switch lower(opt.order)
    case 'lat'
        lonDeg = atan2d(ctx.u(:, 2), ctx.u(:, 1));
        latDeg = asind(ctx.u(:, 3));
        [~, ord] = sortrows([latDeg, lonDeg]);
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
P     = cl(zeros(N, N));       % pre-allocated: no reallocation, no doubling
prev  = 0;
tStart = tic;
nPair  = 0;

for i = 1:nb
    ik = (starts(i):ends(i)).';
    nk = numel(ik);
    Aii = cl(tr_cov_pairs(ctx.Gcov, xyz(ik, :), r(ik), xyz(ik, :), r(ik)));
    nPair = nPair + nk^2;
    dg = (1:(nk+1):nk^2).';
    Aii(dg) = Aii(dg) + ds(ik);

    if i == 1
        S  = (Aii + Aii.') / 2;
        Si = S \ cl(eye(nk, 'like', S));
        alpha(ik)    = Si * ys(ik);
        P(1:nk, 1:nk) = Si;
        prev = nk;
    else
        ip  = (1:prev).';
        ipk = (prev+1:prev+nk).';
        Aij = cl(tr_cov_pairs(ctx.Gcov, xyz(ik, :), r(ik), xyz(ip, :), r(ip)));  % nk x prev
        nPair = nPair + nk * prev;

        v  = ys(ik) - Aij * alpha(ip);
        % PA = P(ip,ip) * Aij' evaluated row block by row block.  Writing
        % P(ip,ip)*Aij' would copy the whole leading block of P, i.e. up to
        % (N - n_b)^2 elements = as much memory as P itself, and the peak
        % memory of the recursion would become 2 N^2 instead of N^2.
        PA   = zeros(prev, nk, 'like', P);
        bwPA = 1024;
        for c0 = 1:bwPA:prev
            cc = (c0:min(c0 + bwPA - 1, prev)).';
            PA(cc, :) = P(ip(cc), ip) * Aij.';
        end

        S  = Aii - Aij * PA;
        S  = (S + S.') / 2;                     % symmetrise, damp rounding
        Si = S \ cl(eye(nk, 'like', S));

        alpha(ip) = alpha(ip) - PA * (Si * v);
        alpha(ik) = Si * v;

        K = PA * Si;
        bw = 2048;                              % row blocks of the gain update
        for c0 = 1:bw:prev
            cc = (c0:min(c0 + bw - 1, prev)).';
            P(ip(cc), ip) = P(ip(cc), ip) + (PA(cc, :) * Si) * PA.';
        end
        P(ip, ipk)  = -K;
        P(ipk, ip)  = -K.';
        P(ipk, ipk) = Si;
        prev = prev + nk;
    end

    if opt.verbose && (mod(i, max(1, ceil(nb / 10))) == 0 || i == nb)
        fprintf('   seq LSC %3d/%3d  t = %6.1f s\n', i, nb, toc(tStart));
    end
end

alpha = alpha(back);
E = -dvar .* double(alpha);
info = struct('time', toc(tStart), 'nBlock', nb, 'blockSize', opt.blockSize, ...
              'nPairEval', nPair, 'pBytes', pBytes, 'single', opt.single, ...
              'exact', true, 'stopReason', 'single-pass-exact');
end
