function D = load_colorado(cfg)
%LOAD_COLORADO  Load the Colorado test field of the 1-cm geoid experiment.
%
%   D = LOAD_COLORADO(cfg) reads the residual gravity disturbances and the
%   GNSS/levelling check points that are distributed with this package
%   (data/colorado) and returns them in the form expected by the solvers.
%
%   The residual field is the one resulting from the remove-compute-restore
%   reduction described in the paper: XGM2019e (degree 719) and ERTM2160 have
%   been removed from the gravity data; the GNSS/levelling height anomalies
%   are reduced in the same way, so that the collocation estimate is a
%   residual that is added back to the reference models before validation.
%
%   Output
%     D.region      'colorado'
%     D.BLh         N x 3  [geodetic latitude, longitude (deg), ellipsoidal height (m)]
%     D.y           N x 1  residual gravity disturbance (mGal)
%     D.n           N
%     D.xyzr        N x 4  [X Y Z r]  geocentric cartesian coordinates and radius
%     D.u           N x 3  unit vector of the observation point
%     D.r           N x 1  geocentric radius (m)
%     D.gnssBLh     222 x 3  GNSS/levelling check points
%     D.gnssObs     222 x 1  height anomaly from GNSS/levelling (m)
%     D.gnssModel   222 x 1  height anomaly of the reference models (m)
%     D.gnssGnorm   222 x 1  normal gravity at the check points (mGal)
%     D.gnssDatum   datum offset of the reference models (m), 0 for this reduction
%
%   The GNSS/levelling points are used for validation only: no parameter of the
%   covariance model, no noise variance and no stopping threshold is derived
%   from them.
%
%   See also SELECT_SUBREGION, GSLSC_CONFIG.

dirData = fullfile(cfg.dataDir, 'colorado');

% ---- gravity observations -------------------------------------------------
f = fullfile(dirData, 'g_res_BLh_xgm2019e.txt');
if ~exist(f, 'file')
    error('load_colorado:missingData', ...
        ['Observation file not found:\n  %s\n' ...
         'The package expects data/colorado/g_res_BLh_xgm2019e.txt, which is ' ...
         'distributed with the repository.'], f);
end
g = load(f);
D.region = 'colorado';
D.BLh    = g(:, 1:3);
D.y      = g(:, 4);
D.n      = size(g, 1);

% ---- GNSS/levelling check points ------------------------------------------
% Preferred source: the XGM2019e-consistent reduction distributed with the
% package.  Fallback: rebuild it from the GSVS17 profiles.
fmat = fullfile(dirData, 'gnss_zeta_xgm2019e.mat');
if exist(fmat, 'file')
    S    = load(fmat);
    D.gnssBLh   = S.BLhG;
    D.gnssObs   = S.zetaObs(:);
    D.gnssModel = S.zModel(:);
    D.gnssDatum = 0;
else
    gs  = load(fullfile(dirData, 'GSVS17.txt'));
    dhm = load(fullfile(dirData, 'DHM_GSVS17.txt'));
    D.gnssBLh   = gs(:, 1:3);
    D.gnssObs   = (gs(:, 3) - gs(:, 4)) - dhm(:);
    D.gnssModel = load(fullfile(dirData, 'gnss_leveling_model_h0_BLh.txt'));
    D.gnssDatum = 0.58;
end

% Normal gravity at the check points.  It is evaluated at the orthometric
% height (the telluroid), which is how the geoid height is obtained from the
% disturbing potential.  GSVS17.txt holds [lat lon h H]; column 4 is H.
gsFile = fullfile(dirData, 'GSVS17.txt');
E  = make_Ellipsoid(cfg.ellipsoid);
if exist(gsFile, 'file')
    gs = load(gsFile);
    if size(gs, 1) ~= size(D.gnssBLh, 1)
        error('load_colorado:gnssMismatch', ...
            'GSVS17.txt has %d points but the height-anomaly file has %d.', ...
            size(gs, 1), size(D.gnssBLh, 1));
    end
    Hortho = gs(:, 4);
else
    Hortho = D.gnssBLh(:, 3);
end
D.gnssGnorm = normal_gravity_confocal(1, E, D.gnssBLh(:, 1), Hortho);

% ---- geocentric coordinates ----------------------------------------------
D.xyzr = BLh_XYZR(D.BLh);
D.u    = D.xyzr(:, 1:3) ./ D.xyzr(:, 4);
D.r    = D.xyzr(:, 4);
end
