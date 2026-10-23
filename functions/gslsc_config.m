function cfg = gslsc_config()
%GSLSC_CONFIG  Paths and default parameters of the GS-LSC package.
%
%   cfg = GSLSC_CONFIG() returns a structure with the directory layout of this
%   repository and with the numerical configuration used in the paper
%
%       Regional gravity field modelling by least-squares collocation with a
%       block Gauss-Seidel iterative solver
%       (Zhang, Chang & Chang, submitted to Computers & Geosciences)
%
%   Every script of the package starts with
%
%       cfg = gslsc_config();
%       addpath(cfg.funDir);
%
%   so that the package can be moved or cloned anywhere: all paths are derived
%   from the location of this file.
%
%   Fields
%     .root       repository root
%     .funDir     functions/        (solver and covariance library)
%     .dataDir    data/            (Colorado test field)
%     .resDir     results/         (created on demand; not version controlled)
%     .figDir     figures/         (created on demand; not version controlled)
%     .cov        Tscherning-Rapp covariance parameters used in the paper
%     .dvar       observation noise variance (mGal^2)
%     .nLat,.nLon latitude/longitude limits of the sub-region of Table 2 (deg)
%     .fun,.crossFun  functional codes of the observation type and of the
%                     observation-to-height-anomaly cross-covariance
%                 ( 11 T/T, 22 gravity anomaly, 33 gravity disturbance,
%                   12 T/gravity anomaly, 13 T/gravity disturbance )
%
%   See also LOAD_COLORADO, SELECT_SUBREGION, LSC_STANDARD,
%   LSC_SEQUENTIAL, LSC_BLOCK_GAUSS_SEIDEL.

cfg.root    = fileparts(fileparts(mfilename('fullpath')));
cfg.funDir  = fullfile(cfg.root, 'functions');
cfg.dataDir = fullfile(cfg.root, 'data');
cfg.resDir  = fullfile(cfg.root, 'results');
cfg.figDir  = fullfile(cfg.root, 'figures');

for d = {cfg.resDir, cfg.figDir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end
addpath(cfg.funDir);

% ---------------------------------------------------------------------
% Covariance model and noise (Table 2 of the paper)
% ---------------------------------------------------------------------
% Tscherning-Rapp kernel, frequency band n > M retained after the removal of
% the reference model XGM2019e (M = 719) and of the topographic model ERTM2160.
cfg.cov = struct('A', 21, 'B', 24, 'M', 719, 'RB', 6369500, 'psiMax', 12);

% Observation noise variance D (mGal^2); the paper uses D = 16 mGal^2, i.e. a
% standard deviation of 4 mGal, consistent with the airborne/terrestrial mix of
% the 1-cm geoid experiment.  Section 4.2 of the paper repeats the experiment
% with cfg.dvar = 1 mGal^2 to illustrate the dependence on the condition number.
cfg.dvar = 16;

% ---------------------------------------------------------------------
% Test field: Colorado sub-region (Table 2)
% ---------------------------------------------------------------------
% The box below is the sub-region used in the paper.  Note that the box
% contains 50 902 of the 93 588 observations of the full Colorado field; the
% paper rounds this to N = 5 x 10^4.  Set cfg.nMax = 50000 to draw a fixed-seed
% subsample of that size instead (see SELECT_SUBREGION).
cfg.nLat = [34.972 39.171];
cfg.nLon = [-108.203 -104.003];
cfg.nMax = [];

% residual gravity disturbance (T - gravity disturbance cross-covariance is 13)
cfg.fun      = 33;
cfg.crossFun = 13;

% Earth model used by the normal gravity field
cfg.ellipsoid = 'GRS80';

% ---------------------------------------------------------------------
% Solver defaults
% ---------------------------------------------------------------------
cfg.blockSize = 2000;      % nb of the block solvers (observations per block)
cfg.precision = 'single';  % precision of the standard LSC matrix (9.3 GB at N=5e4)
end
