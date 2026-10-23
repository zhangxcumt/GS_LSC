function p = legendreX(t, Nmax)
%LEGENDREX  Legendre polynomials P_0(t) ... P_Nmax(t).
%
%   p = LEGENDREX(t, Nmax)
%
%   Input
%     t    - argument, any shape (used as a row of values)
%     Nmax - highest degree
%   Output
%     p    - (Nmax+1) x numel(t); row k holds P_{k-1}(t)
%
%   Standard three-term recurrence
%   (n+1) P_{n+1} = (2n+1) t P_n - n P_{n-1}.

t = t(:).';
p = zeros(Nmax + 1, numel(t));
p(1, :) = 1;
if Nmax >= 1, p(2, :) = t; end
if Nmax >= 2, p(3, :) = 0.5 * (3 * t.^2 - 1); end
if Nmax >= 3, p(4, :) = 0.5 * (5 * t.^3 - 3 * t); end
% Three-term recurrence: n P_n = (2n-1) t P_{n-1} - (n-1) P_{n-2};
% row k of p holds P_{k-1}.
for n = 4:Nmax
    p(n + 1, :) = ((2 * n - 1) * t .* p(n, :) - (n - 1) * p(n - 1, :)) / n;
end
end
