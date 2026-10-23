function [FB] = cal_FB(s, t, B)
%CAL_FB  Closed form of the Tscherning-Rapp kernel F_B(s,t).
%
%   F_B(s,t) = sum_{n=3}^{inf} s^(n+1) * P_n(t) / (n + B)
%
%   s = RB^2/(r_P r_Q) is the Bjerhammar radius ratio, t = cos(psi) the cosine
%   of the spherical distance and P_n the (unnormalised) Legendre polynomials.
%   The terms n = 0,1,2 are subtracted inside the function, so the return value
%   is exactly the infinite series starting at n = 3 (Tscherning & Rapp, 1974).
%
%   Input
%     s - radius ratio (scalar or array)
%     t - cosine of the spherical distance (same shape as s, or broadcastable)
%     B - model parameter B (-2, -1, 0, or a positive integer; B = 24 here)
%
%   Output
%     FB - F_B(s,t), same shape as the input
%
%   Reference: Tscherning C.C., Rapp R.H. (1974) Closed covariance expressions
%   for gravity anomalies, geoid undulations and deflections of the vertical
%   implied by anomaly degree variance models. OSU Report 208.

L = sqrt(1 - 2 .* s .* t + s.^2);
N = 1 + L - s .* t;
M = 1 - L - s .* t;

% low-degree terms (n = 0,1,2): p_n = s^(n+1) P_n(t)
p0 = s;
p1 = s.^2 .* t;
p2 = s.^3 .* (1.5 .* t.^2 - 0.5);

switch B
    case -2
        lowdeg = 0;
        FB = 0.5 .* s .* (1 + 3 .* s .* t) .* M + p2 .* log(2 ./ N) ...
             + 0.25 .* s.^3 .* (1 - t.^2);
    case -1
        lowdeg = p2;
        FB = s .* (M + s .* t .* log(2 ./ N));
    case 0
        lowdeg = p1 + p2 / 2;
        FB = s .* log(2 ./ N);
    case 1
        lowdeg = p0 / 1 + p1 / 2 + p2 / 3;
        FB = log(1 + 2 .* s ./ (1 - s + L));
    case 2
        lowdeg = p0 / 2 + p1 / 3 + p2 / 4;
        F1 = log(1 + 2 .* s ./ (1 - s + L));
        FB = (L - 1 + t .* F1) ./ s;
    otherwise
        % recurrence: from F_1 and F_2 up to the required order
        lowdeg = p0 / B + p1 / (B + 1) + p2 / (B + 2);
        F1 = log(1 + 2 .* s ./ (1 - s + L));      % F_1
        F2 = (L - 1 + t .* F1) ./ s;              % F_2
        if B >= 3
            FB = F2;                              % B = 3: the loop runs once
            for BB = 3:B
                FB = (L + (2 * BB - 3) .* t .* F2 - (BB - 2) .* F1 ./ s) ./ (s * (BB - 1));
                F1 = F2;
                F2 = FB;
            end
        else
            error('cal_FB: B must be an integer not smaller than -2.');
        end
end

FB = FB - lowdeg;
end
