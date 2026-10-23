function G = tr_cov_grid(P, functional, opt)
%TR_COV_GRID  Tabulate the Tscherning-Rapp covariance on the (t, s) plane.
%
%   The covariance function depends only on the radius ratio
%   s = RB^2/(r_P r_Q) and on the cosine of the spherical distance t = cos(psi).
%   Tabulating it once on the (t, s) plane reduces every later point pair to a
%   two-dimensional interpolation.  This is what makes N ~ 10^5 observations
%   tractable: the closed form costs a Legendre series of M terms, the table
%   costs a bilinear interpolation.
%
%   The s grid is uniform in u = 1 - s on a logarithmic scale.  For the band
%   n > M the covariance behaves as C ~ sum s^n c_n with an effective degree
%   n_eff ~ 1/u, so its relative variation scale is ~ u and a logarithmic grid
%   reaches the same accuracy with far fewer nodes than a uniform one.
%
%   Input
%     P          - parameter struct A, B, M, RB
%     functional - functional code, see TR_COV_CLOSED
%     opt        - .psiMax (12)  maximum spherical distance (deg)
%                  .resT (max(4000,10*M))  number of t nodes
%                  .resS (1201)            number of s nodes
%                  .rRef   radii of the observation points, for the s range
%                  .sRange explicit [sLo sHi] instead of the adaptive range
%                  .sMargin (0.05) relative margin added to the s range
%
%   Output
%     G - .F (griddedInterpolant, call F(tq,sq)), .tvec, .svec, .info

if nargin < 3, opt = struct(); end
if ~isfield(opt,'psiMax') || isempty(opt.psiMax), opt.psiMax = 12; end
if ~isfield(opt,'resT')   || isempty(opt.resT),   opt.resT   = max(4000, 10*P.M); end
if ~isfield(opt,'resS')   || isempty(opt.resS),   opt.resS   = 1201; end
if ~isfield(opt,'sMargin')|| isempty(opt.sMargin),opt.sMargin= 0.05; end

RB = P.RB; A = P.A; B = P.B; M = P.M;

% ---- range of s, parametrised by u = 1 - s ----
if isfield(opt,'sRange') && ~isempty(opt.sRange)
    sLo = opt.sRange(1); sHi = opt.sRange(2);
else
    if isfield(opt,'rRef') && ~isempty(opt.rRef)
        rr = opt.rRef(:);
        sLo = (RB / max(rr))^2;
        sHi = (RB / min(rr))^2;
    else
        sLo = 0.996; sHi = 1.0;
    end
    span = max(sHi - sLo, 1e-9);
    sLo = max(1e-9, sLo - opt.sMargin * span);
    sHi = sHi + opt.sMargin * span;
end
uLo = max(1 - sHi, 1e-12);
uHi = 1 - sLo;
uvec = logspace(log10(uHi), log10(uLo), opt.resS);   % decreasing
svec = 1 - uvec;                                      % increasing
svec(end) = min(svec(end), 1 - 1e-12);

% ---- t grid: from cos(psiMax) up to 1 ----
psi = linspace(0, opt.psiMax, opt.resT);
tvec = fliplr(cosd(psi));
tvec(1) = max(tvec(1), -1);

% ---- Legendre basis P_3..P_M ----
Pn = legendreX(tvec, M);
Pn = Pn(4:end, :);

nn = (3:M).';
switch functional
    case 11
        coef = A * RB^2 ./ ((nn - 1) .* (nn - 2) .* (nn + B));        pows = nn + 1;
    case 22
        coef = A * (nn - 1) ./ ((nn - 2) .* (nn + B));                pows = nn + 2;
    case 33
        coef = A * (nn + 1).^2 ./ ((nn - 1) .* (nn - 2) .* (nn + B)); pows = nn + 2;
    case 12
        coef = A ./ ((nn - 2) .* (nn + B));                           pows = nn + 2;
    case 13
        coef = A * (nn + 1) ./ ((nn - 1) .* (nn - 2) .* (nn + B));    pows = nn + 2;
    otherwise
        error('tr_cov_grid: invalid functional code.');
end

resS = numel(svec); resT = numel(tvec);
Cmat = zeros(resS, resT);
for i = 1:resS
    Cmat(i, :) = tr_cov_closed(svec(i), tvec, P, functional, Pn);
end

G = struct();
G.F = griddedInterpolant({tvec(:), svec(:)}, Cmat.', 'linear', 'linear');
G.tvec = tvec; G.svec = svec;
G.info = struct('functional', functional, 'A', A, 'B', B, 'M', M, 'RB', RB, ...
                'sRange', [svec(1) svec(end)], 'uRange', [min(uvec) max(uvec)], ...
                'resT', resT, 'resS', resS, 'psiMax', opt.psiMax);
end
