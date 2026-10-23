function sel = select_subregion(D, cfg, nMax, seed)
%SELECT_SUBREGION  Indices of the Colorado sub-region used in the paper.
%
%   sel = SELECT_SUBREGION(D, cfg) returns the indices of the observations that
%   lie inside the box of Table 2 of the paper,
%
%       latitude   34.972 deg N  ...  39.171 deg N
%       longitude 108.203 deg W  ... 104.003 deg W
%
%   The box contains 50 902 of the 93 588 observations.  The paper quotes
%   N = 5 x 10^4; to reproduce exactly that count use
%
%       sel = SELECT_SUBREGION(D, cfg, 50000, 20260923)
%
%   which draws nMax observations from the box with a fixed seed (so the
%   selection is reproducible and identical on every run and every platform).
%
%   Unlike spatial decimation (e.g. keeping every k-th point), the sub-region
%   retains the full observation density inside the area: the point spacing is
%   the one of the test field itself.
%
%   Input
%     D    - structure returned by LOAD_COLORADO
%     cfg  - structure returned by GSLSC_CONFIG (uses .nLat, .nLon)
%     nMax - (optional) number of observations to keep; [] keeps the whole box
%     seed - (optional) seed of the random draw, default 20260923
%
%   Output
%     sel  - column vector of indices into D.y, sorted ascending
%
%   See also LOAD_COLORADO.

if nargin < 3 || isempty(nMax), nMax = cfg.nMax; end
if nargin < 4 || isempty(seed), seed = 20260923; end

lat = D.BLh(:, 1);
lon = D.BLh(:, 2);
inside = lat >= cfg.nLat(1) & lat <= cfg.nLat(2) & ...
         lon >= cfg.nLon(1) & lon <= cfg.nLon(2);
sel = find(inside);

if ~isempty(nMax) && nMax < numel(sel)
    rng(seed, 'twister');
    sel = sort(sel(randperm(numel(sel), nMax)));
end
sel = sel(:);
end
