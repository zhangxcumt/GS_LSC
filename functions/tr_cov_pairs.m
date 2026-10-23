function C = tr_cov_pairs(G, xyzP, rP, xyzQ, rQ)
%TR_COV_PAIRS  Cross-covariance matrix from the tabulated covariance function.
%
%   C(k,l) = C_{fg}(P_k, Q_l),  f = functional of the prediction point,
%                               g = functional of the observation point
%
%   Input
%     G    - output of TR_COV_GRID
%     xyzP - nP x 3 geocentric cartesian coordinates of the prediction points
%     rP   - nP x 1 geocentric radii of the prediction points
%     xyzQ - nQ x 3 geocentric cartesian coordinates of the observation points
%     rQ   - nQ x 1 geocentric radii of the observation points
%
%   Output
%     C    - nP x nQ
%
%   For functionals 12 and 13 (cross-covariance of T with gravity) the radius
%   of the prediction point is applied here; TR_COV_GRID.info.functional decides
%   whether this factor is needed.

RB = G.info.RB;
fun = G.info.functional;
if fun == 12 || fun == 13
    scaleP = rP(:);
else
    scaleP = ones(numel(rP), 1);
end

uP = xyzP ./ rP(:);
uQ = xyzQ ./ rQ(:);

t = uP * uQ.';                       % nP x nQ, cos(psi)
s = (RB^2) ./ (rP(:) * rQ(:).');     % nP x nQ

% keep the arguments inside the tabulated ranges (rounding can leave them out)
t = min(max(t, G.tvec(1)), G.tvec(end));
s = min(max(s, G.svec(1)), G.svec(end));

C = G.F(t, s);
C = C .* scaleP;
end
