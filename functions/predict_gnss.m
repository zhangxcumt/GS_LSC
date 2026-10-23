function [zetaRes, S] = predict_gnss(Gcross, coef, ctx, D)
%PREDICT_GNSS  Predict the residual height anomaly at the GNSS/levelling points.
%
%   [zetaRes, S] = PREDICT_GNSS(Gcross, coef, ctx, D)
%
%   Adds the collocation estimate of the residual field to the height anomaly
%   of the reference models (XGM2019e + ERTM2160) and compares it with the
%   GNSS/levelling height anomalies.
%
%   The mean difference (a datum offset) is REMOVED before the statistics are
%   computed and reported separately in S.datum: the collocation estimate is
%   not allowed to absorb it, and no parameter of the model is derived from the
%   check points.
%
%   Input
%     Gcross - T - gravity cross-covariance table (ctx.Gcross)
%     coef   - N x 1 collocation coefficients
%     ctx    - context of the observations
%     D      - structure from LOAD_COLORADO
%
%   Output
%     zetaRes - n x 1 residual height anomaly at the check points (m)
%     S       - .n .std .rms .max .min .datum .err
%               std, rms, max, min and err are in metres, after removing the
%               datum offset; err is the vector of errors used by the figures.
%
%   See also SIGNAL_PREDICT, LOAD_COLORADO.

xyzG = BLh_XYZR(D.gnssBLh);
zetaRes = signal_predict(Gcross, coef, ctx, xyzG(:, 1:3), xyzG(:, 4), ...
    struct('gnorm', D.gnssGnorm));

full = zetaRes + D.gnssModel + D.gnssDatum;
gd   = full - D.gnssObs;
datum = mean(gd);          % datum offset of the reference models, reported only
gd    = gd - datum;

S = struct('n', numel(gd), 'std', std(gd), 'rms', sqrt(mean(gd.^2)), ...
           'max', max(gd), 'min', min(gd), 'datum', datum, 'err', gd);
end
