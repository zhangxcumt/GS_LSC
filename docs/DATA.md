# Data distributed with the package

Directory: `data/colorado/`

| file | size | content |
|---|---:|---|
| `g_res_BLh_xgm2019e.txt` | 4.1 MB | 93 588 residual gravity disturbances: `lat lon h y` (deg, deg, m, mGal) |
| `GSVS17.txt` | 18 kB | 222 GNSS/levelling points of the GSVS17 line: `lat lon h H` (deg, deg, m, m) |
| `DHM_GSVS17.txt` | 4 kB | digital height model correction along GSVS17 (m), one value per point |
| `gnss_leveling_model_h0_BLh.txt` | 4 kB | height anomaly of the reference models at the check points (m) |
| `gnss_zeta_xgm2019e.mat` | 14 kB | the same check points reduced consistently with the gravity data: `BLhG` (222 × 3), `zetaObs`, `zModel` (m) |

`load_colorado.m` reads `g_res_BLh_xgm2019e.txt` and, for the check points,
`gnss_zeta_xgm2019e.mat` when it exists and the three ASCII files otherwise.

## Provenance

The Colorado test field is the public test field of the **1-cm geoid
experiment**:

* Wang Y.M., Sánchez L., Ågren J., Huang J., Forsberg R. (2021) Colorado geoid
  computation experiment: overview and summary. *Journal of Geodesy* 95, 127.
* The GSVS17 GNSS/levelling profile is the official validation data set of that
  experiment.

The residual field distributed here was obtained within the
remove-compute-restore framework with

* **XGM2019e** as the reference global gravity field model, removed to degree
  719 — Zingerle P., Pail R., Gruber T., Oikonomidou X. (2020) The combined
  global gravity field model XGM2019e. *Journal of Geodesy* 94;
* **ERTM2160** as the topographic gravity field model — Hirt C., Kuhn M.,
  Claessens S., Pail R., Seitz K., Gruber T. (2014) Study of the Earth's
  short-scale gravity field using the ERTM2160 gravity model.
  *Computers & Geosciences* 73, 71–80.

The residual gravity disturbances are the terrestrial and airborne observations
of the experiment, after the reduction; the airborne data are the
Gauss-filtered, decimated data of the experiment protocol.

## Terms of use

The data originate from the 1-cm geoid experiment and are redistributed here
**only to make the computations of the paper reproducible**. Users must respect
the terms of the original providers. If you use the data in a publication,
cite Wang et al. (2021) and the data providers of the experiment.

The code in this repository is licensed under the MIT license (see `LICENSE`);
that license covers the code, not the data.

## Validation protocol

The 222 GNSS/levelling points are used **for validation only**:

* no covariance parameter is estimated from them (`gslsc_config.m` fixes
  A, B, M and R_B from the gravity data alone),
* the noise variance `D` is a prescribed constant, not tuned on the check
  points,
* the stopping threshold of the iterative solver refers to the change of the
  predicted quantity, never to its agreement with the check points,
* the mean difference between model and check points (the datum offset, about
  0.87 m for this reduction) is removed before the statistics are computed and
  reported separately, e.g. in `predict_gnss.m` and in Table 3 of the paper.

## Regenerating the data

Recomputing the residual field from the original gravity data requires the
XGM2019e and ERTM2160 coefficients, which are not redistributed here because of
their size. The scripts are not part of this package; the residual field is
distributed as a ready-to-use ASCII file, which is what the solvers need.
