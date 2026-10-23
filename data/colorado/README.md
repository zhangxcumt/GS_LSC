# Colorado test field

This directory holds the test field used in the paper. See `docs/DATA.md` at the
root of the repository for the provenance, the terms of use and the validation
protocol.

| file | content |
|---|---|
| `g_res_BLh_xgm2019e.txt` | 93 588 residual gravity disturbances `lat lon h y` (deg, deg, m, mGal), XGM2019e (719) + ERTM2160 removed |
| `GSVS17.txt` | 222 GNSS/levelling points `lat lon h H` (deg, deg, m, m) |
| `DHM_GSVS17.txt` | digital height model correction along GSVS17 (m) |
| `gnss_leveling_model_h0_BLh.txt` | height anomaly of the reference models at the check points (m) |
| `gnss_zeta_xgm2019e.mat` | check points reduced consistently with the gravity data: `BLhG`, `zetaObs`, `zModel` (m) |

Source: the 1-cm geoid experiment (Wang et al., 2021, *J. Geod.* 95, 127).
Redistributed for reproducibility only; cite the original publication.
