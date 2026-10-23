function C = tr_cov_closed(s, t, P, functional, Pn)
%TR_COV_CLOSED  Band-limited closed expression of the Tscherning-Rapp covariance.
%
%   C = A * sum_{n=M+1}^{inf} s^p(n) * c_n * P_n(t)
%
%   Evaluated as "infinite closed form minus finite series": the closed form of
%   CAL_FB gives the complete sum n = 3..inf, from which the explicit terms
%   n = 3..M are subtracted.  The result is therefore exactly the covariance of
%   the frequency band n > M, i.e. of the residual field left after the
%   reference model has been removed to degree M.
%
%   Input
%     s          - Bjerhammar radius ratio s = RB^2/(r_P r_Q), scalar or vector
%     t          - cosine of the spherical distance, row or column vector
%     P          - parameter struct with fields A, B, M, RB
%     functional - functional combination:
%                    11  T - T                22  gravity anomaly - gravity anomaly
%                    33  gravity disturbance - gravity disturbance
%                    12  T - gravity anomaly  13  T - gravity disturbance
%                  Functionals 12 and 13 return the cross-covariance between T
%                  and gravity; it must be multiplied by the radius of the
%                  prediction point, which TR_COV_PAIRS does automatically.
%     Pn         - (optional) precomputed Legendre matrix P_3..P_M, size
%                  (M-2) x numel(t)
%
%   Coefficient c_n and power p(n) of each functional:
%     11 : c_n = A RB^2 /((n-1)(n-2)(n+B)),      p = n+1
%     22 : c_n = A (n-1)/((n-2)(n+B)),           p = n+2
%     33 : c_n = A (n+1)^2/((n-1)(n-2)(n+B)),    p = n+2
%     12 : c_n = A /((n-2)(n+B)),                p = n+2
%     13 : c_n = A (n+1)/((n-1)(n-2)(n+B)),      p = n+2
%
%   Reference: Tscherning & Rapp (1974).

A = P.A; B = P.B; M = P.M; RB = P.RB;
t = t(:).';
scalarS = isscalar(s);
if scalarS
    sv = s;
else
    sv = reshape(s, 1, []);
    if numel(sv) ~= numel(t)
        error('tr_cov_closed: s and t must have the same length (or s scalar).');
    end
end
if (nargin < 5) || isempty(Pn)
    Pn = legendreX(t, M);
    Pn = Pn(4:end, :);
end
nn = (3:M).';

switch functional
    case 11
        cf = A * RB^2 / ((B + 1) * (B + 2)) * ...
             ((B + 1) * cal_FB(sv, t, -2) - (B + 2) * cal_FB(sv, t, -1) + cal_FB(sv, t, B));
        coef = A * RB^2 ./ ((nn - 1) .* (nn - 2) .* (nn + B));  pows = nn + 1;
    case 22
        cf = A * sv / (B + 2) .* (cal_FB(sv, t, -2) + (B + 1) * cal_FB(sv, t, B));
        coef = A * (nn - 1) ./ ((nn - 2) .* (nn + B));          pows = nn + 2;
    case 33
        cf = A * sv / ((B + 1) * (B + 2)) .* (9 * (B + 1) * cal_FB(sv, t, -2) ...
             - 4 * (B + 2) * cal_FB(sv, t, -1) + (B - 1)^2 * cal_FB(sv, t, B));
        coef = A * (nn + 1).^2 ./ ((nn - 1) .* (nn - 2) .* (nn + B)); pows = nn + 2;
    case 12
        cf = A * sv / (B + 2) .* (cal_FB(sv, t, -2) - cal_FB(sv, t, B));
        coef = A ./ ((nn - 2) .* (nn + B));                     pows = nn + 2;
    case 13
        cf = A * sv / ((B + 1) * (B + 2)) .* (-2 * (B + 2) * cal_FB(sv, t, -1) ...
             + 3 * (B + 1) * cal_FB(sv, t, -2) - (B - 1) * cal_FB(sv, t, B));
        coef = A * (nn + 1) ./ ((nn - 1) .* (nn - 2) .* (nn + B));    pows = nn + 2;
    otherwise
        error('tr_cov_closed: functional must be one of 11, 22, 33, 12, 13.');
end

if scalarS
    band = (coef .* sv.^pows).' * Pn;
else
    band = sum((coef .* sv.^pows) .* Pn, 1);
end
C = cf - band;
end
