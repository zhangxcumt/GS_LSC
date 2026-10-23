# GS-LSC — least-squares collocation with a block Gauss-Seidel iterative solver

MATLAB code accompanying

> **Regional gravity field modelling by least-squares collocation with a block
> Gauss-Seidel iterative solver**
> Xun Zhang, Lubin Chang\*, Guobin Chang
> *Computers & Geosciences* (submitted)

Corresponding author: Lubin Chang <1410021052@nue.edu.cn>

The package contains the three collocation solvers that are compared in the
paper, the Tscherning–Rapp covariance model they use, the Colorado test field of
the 1-cm geoid experiment, and the scripts that reproduce Tables 3–4 and
Figs 1, 4 and 5.

---

## 1. What the code does

Least-squares collocation solves

```
(C + D) alpha = y ,        s(P) = sum_j alpha_j C_sg(P, Q_j)
```

with `C` the N × N covariance matrix of the observations, `D` the observation
noise covariance, and `alpha` the vector of collocation coefficients. The three
solvers differ only in how they solve that system:

| file | method | matrix inverted | working storage |
|---|---|---|---|
| `functions/lsc_standard.m` | standard LSC, Cholesky factorisation of `A = C + D` | N × N | O(N²) |
| `functions/lsc_sequential.m` | sequential (stepwise) LSC, Moritz (1973) / Tscherning (1974) | n_b × n_b | O(N²) — the recursion state `P` |
| `functions/lsc_block_gauss_seidel.m` | block Gauss-Seidel iteration on `A alpha = y` | n_b × n_b | O(n_b² + N) |

The point of the paper is that sequential collocation reduces the size of the
matrix **to be inverted**, not the total storage, and that a block Gauss-Seidel
iteration is the only one of the three whose storage is independent of `N`,
at the price of re-evaluating all N² covariance elements at every sweep.

Two further solvers are provided as reference points: `functions/lsc_matvec.m`
(matrix–vector product with `A` without forming `A`) and
`functions/lsc_pcg.m` (Jacobi-preconditioned conjugate gradients).

---

## 2. Requirements

* MATLAB **R2018b or newer**. No toolbox is required.
  (The paper was produced with R2026a; `griddedInterpolant`, `chol` on single
  precision and `sortrows` are the only non-trivial dependencies.)
* About **20 GB of RAM** for the `standard_lsc` experiment at N = 50 902, and
  **128 MB** for the block Gauss-Seidel experiment.
* Any operating system; paths are derived from the location of this file.

---

## 3. Quick test

The quick test is the one referred to in the *Computer Code Availability*
section of the manuscript. It takes about **one minute** and needs no
additional data — the Colorado test field is distributed in `data/`.

From MATLAB, in the repository root:

```matlab
cd examples
run_quick_test
```

or, without changing directory:

```matlab
run('examples/run_quick_test.m')
```

It runs two examples:

1. `example_known_solution.m` (a few seconds, no data). A synthetic observation
   vector is built from a known coefficient vector, `y = (C + D) alpha_true`,
   and the three solvers must recover `alpha_true`. Prints the relative errors
   and the tolerances, and stops with an error if any check fails.
2. `example_colorado_subset.m` (about 15 s). The three solvers are run on a
   compact subsample of 5 895 real Colorado observations and the height
   anomalies predicted at the 222 GNSS/levelling check points are compared.

Expected output (values vary slightly with the MATLAB version and the machine):

```
   PASS: the three solvers reproduce the known solution.

   solver                         STD (cm)     time (s)
   standard LSC                   4.0837       0.54
   sequential LSC                 4.0837       2.02
   block Gauss-Seidel LSC         4.0848       13.52

   max |dzeta|, sequential vs standard      : 1.127e-11 cm
   max |dzeta|, block GS  vs standard      : 3.969e-02 cm
   |STD(block GS) - STD(standard)|         : 1.104e-03 cm

   PASS: the three solvers predict the same height anomaly.

################ QUICK TEST PASSED ################
```

The `0.04 cm` of the block Gauss-Seidel row is not an error of the solver: it is
the last change of the predicted quantity, accumulated over the final sweeps
before the 0.01 cm stopping threshold was reached. See Section 4.1 of the paper.

---

## 4. Reproducing the paper

All experiment scripts are in `scripts/` and start with their own `addpath`, so
they can be launched from anywhere. Results go to `results/`, figures to
`figures/` (both are created on demand and are not version controlled).

| Paper item | Script | Output | Cost at N = 50 902 |
|---|---|---|---|
| Table 3, standard LSC | `run_standard_lsc` | `results/standard_lsc.mat` | ~2 min, 20 GB peak |
| Table 3, sequential LSC | `run_sequential_lsc` | `results/sequential_lsc.mat` | ~25 min, 20 GB |
| Table 3, block Gauss-Seidel | `run_block_gs_lsc` | `results/block_gs_lsc_<nb>.mat` | ~35 min, 128 MB |
| Table 3, assembly | `make_table3` | screen + `results/table3.csv` | s |
| Table 4, sweep cost | `run_table4_sweep_cost` | screen | ~1 min |
| Fig. 1, storage vs N | `make_fig1_storage` | `figures/fig1_storage.*` | s |
| Fig. 4, GSVS17 errors | `make_fig4_gsvs17` | `figures/fig4_gsvs17.*` | s |
| Fig. 5, convergence | `make_fig5_convergence` | `figures/fig5_convergence.*` | ~1.5 h |

Everything at once:

```matlab
cd scripts
run_all
```

To check the whole chain on a smaller sample before committing to the full runs,
set the environment variable `GSLSC_NMAX`; every experiment then uses that many
observations instead of the 50 902 of the sub-region:

```powershell
$env:GSLSC_NMAX = 8000      # PowerShell
```
```
set GSLSC_NMAX=8000         # cmd
export GSLSC_NMAX=8000      # bash
```

A second environment variable, `GSLSC_SEQ_SINGLE=1`, makes `run_sequential_lsc`
keep its recursion state `P` in single precision (4 N² instead of 8 N² bytes).
That is how the sequential solution at N = 50 902 can be produced on a machine
where the 19.3 GiB double-precision state does not fit; the script says so when
it is used.

The numerical configuration used in the paper is the set of defaults in
`functions/gslsc_config.m`:

| quantity | value |
|---|---|
| covariance parameters (A, B, M, R_B) | (21, 24, 719, 6 369 500 m) |
| observation noise variance D | 16 mGal² |
| sub-region | 34.972–39.171 °N, 108.203–104.003 °W |
| observations in the sub-region | 50 902 (the paper quotes 5 × 10⁴) |
| GNSS/levelling check points | 222 |
| block sizes n_b | 2000 and 4000 |
| stopping threshold (predicted quantity) | 0.01 cm |

> **Note on N.** The bounding box printed in Table 2 of the paper contains
> **50 902** observations of the 93 588 of the full Colorado field, and that is
> what the scripts use. To reproduce exactly `N = 50 000`, call
> `select_subregion(D, cfg, 50000)`; the setting `cfg.nMax = 50000` in
> `gslsc_config.m` does the same.

### What the package produces

Running `scripts/run_all.m` on the machine of the paper (Intel Core Ultra 7,
31.4 GB RAM, MATLAB R2026a) with all 50 902 observations gives

| quantity | this package (N = 50 902) | paper (N = 5 × 10⁴) |
|---|---:|---:|
| GNSS/levelling RMS, standard LSC | 2.5177 cm | 2.5690 cm |
| GNSS/levelling STD, standard LSC | 2.5234 cm | 2.5132 cm |
| GNSS/levelling RMS / STD, sequential LSC | 2.5177 / 2.5234 cm | 2.5690 / 2.5132 cm |
| GNSS/levelling STD, block GS n_b = 2000 | 2.5310 cm | 2.5242 cm |
| GNSS/levelling STD, block GS n_b = 4000 | 2.5304 cm | 2.5112 cm |
| self-fit STD, standard LSC | 2.8104 mGal | 2.8310 mGal |
| sweeps, block GS n_b = 2000 | 35 | 40 |
| sweeps, block GS n_b = 4000 | 26 | 33 |
| time, standard LSC | 61.6 s | 53.3 s |
| time, sequential LSC | 298 s (P in single) | 1370.6 s (P in double) |
| time, block GS n_b = 2000 | 1173 s | 1106 s |
| time, block GS n_b = 4000 | 805 s | 865 s |
| relative residual at stop | 1.0 × 10⁻² / 1.1 × 10⁻² | 8.3 × 10⁻³ – 2.0 × 10⁻² |
| storage, standard LSC | 9.65 GB | 9.31 GB |
| storage, block GS | 31 MB / 122 MB | 32 MB / 128 MB |
| **one matrix–vector product with A cached** | **0.167 s** | **0.167 s** |
| **one sweep with the covariance recomputed** | **26.4 s** | **28.49 s** |
| **ratio of the two (Table 4)** | **158** | **171** |

The last three rows reproduce Table 4 of the paper to within 8 %, and the sweep
cost agrees to 7 %. Standard and sequential LSC give **identical** values for
every indicator, as the paper states. The differences of 0.01–0.05 cm in the
block Gauss-Seidel rows are those of a slightly different observation subset
(50 902 instead of 50 000), not of a different solver.

The sequential run was made with `GSLSC_SEQ_SINGLE=1`, i.e. with the recursion
state `P` in single precision (9.65 GB instead of 19.3 GB), because the
double-precision state does not fit next to the operating system on a 32 GB
machine that is in normal use. The coefficients are then accurate to about 1e-6
relative and the statistics above are unchanged; run `run_sequential_lsc`
without that variable on a machine with ~26 GB free to reproduce the
double-precision timings of the paper.

---

## 5. Repository layout

```
functions/          the library
  gslsc_config.m                paths and default parameters
  load_colorado.m               test field and check points
  select_subregion.m            the sub-region of Table 2
  make_context.m                covariance tables for a set of points
  lsc_standard.m                standard LSC (Cholesky)
  lsc_sequential.m              sequential LSC (Moritz/Tscherning recursion)
  lsc_block_gauss_seidel.m      block Gauss-Seidel LSC
  lsc_pcg.m                     PCG reference solver
  lsc_matvec.m                  (C + D) * alpha without forming C
  build_covariance_matrix.m     blocked assembly of A = C + D
  signal_predict.m              prediction from the coefficients
  predict_gnss.m                validation at the GNSS/levelling points
  storage_model.m               analytic storage of the three solvers
  setup_experiment.m            common input of the experiment scripts
  plot_fig5.m                   convergence figure from the recorded history
  tr_cov_closed.m               Tscherning-Rapp covariance, closed form
  tr_cov_grid.m                 tabulation on the (t, s) plane
  tr_cov_pairs.m                covariance of point pairs from the table
  cal_FB.m, legendreX.m         closed-form kernel and Legendre polynomials
  normal_gravity_confocal.m, make_Ellipsoid.m, BLh_XYZR.m
examples/           quick test and examples (see Section 3)
scripts/            experiments producing Tables 3-4 and Figs 1, 4, 5
data/colorado/      test field and GNSS/levelling check points
docs/               algorithm notes, data provenance, availability statement
results/, figures/  generated; not version controlled
```

---

## 6. Data

`data/colorado` holds the Colorado test field of the 1-cm geoid experiment
(Wang et al., 2021, *J. Geod.* 95, 127), reduced with XGM2019e (degree 719) and
ERTM2160. The GNSS/levelling points of the GSVS17 line are used **for validation
only**: no covariance parameter, no noise variance and no stopping threshold is
derived from them. See `docs/DATA.md` for the file-by-file description and for
the terms of use.

---

## 7. Numerical notes

* `lsc_sequential` works in the dual and keeps the recursion state `P` of order
  N × N; the function refuses to run when `P` would exceed `opt.maxMemGB`
  (26 GB by default), because this is exactly the limitation the paper
  quantifies. The gain update is accumulated in row blocks, using
  `K S K' = (P A) S^{-1} (P A)'`, so that no `prev × prev` temporary is created.
* `lsc_block_gauss_seidel` carries no state across sweeps: the only unknown is
  `alpha`, and the covariance blocks are recomputed at every sweep. Its stopping
  rule is the change of the **predicted** quantity, not of the residual; pass a
  handle through `opt.thetaFun` (see `run_block_gs_lsc.m`). The residual at the
  stopping point is reported in `info.relResid` for information only.
* The covariance is never evaluated in closed form inside the solvers: it is
  tabulated once on the (t, s) plane by `tr_cov_grid` and interpolated by
  `tr_cov_pairs`, which is what makes O(N²) covariance work affordable.
* All random splits use fixed seeds, so the results are reproducible.

---

## 8. License and citation

MIT license, see `LICENSE`. If you use this code, please cite the article above
and, once it is archived, the DOI of this repository (see `CITATION.cff`).

Contact: Lubin Chang <1410021052@nue.edu.cn>, Xun Zhang <zhangxun@cumt.edu.cn>.
