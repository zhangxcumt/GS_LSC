# Algorithms

This note collects the equations implemented by the three solvers, the
complexity of each, and the implementation choices that matter for
reproducing the numbers of the paper. Notation follows Section 2 of the paper.

## 0. The collocation problem

Given N observations `y` with covariance matrix `C` and noise covariance `D`,

```
A alpha = y ,      A = C + D ,
s(P) = sum_j alpha_j C_sg(P, Q_j) ,      zeta(P) = s(P) / gamma(P)
```

`s` is the residual disturbing potential, `zeta` the residual height anomaly and
`gamma` the normal gravity at the prediction point. The observations are residual
gravity disturbances (mGal) after the removal of XGM2019e (degree 719) and
ERTM2160.

The covariance is the band-limited Tscherning–Rapp model

```
C_fg(P,Q) = sum_{n=M+1}^{inf} A  c_n(f,g)  s^(p_n)  P_n(cos psi) ,
s = R_B^2 / (r_P r_Q) ,
```

evaluated by `tr_cov_closed.m` as *closed form minus the finite series
`n = 3..M`*, and tabulated on the `(t = cos psi, s)` plane by `tr_cov_grid.m`.
`tr_cov_pairs.m` then returns any covariance block by bilinear interpolation.

## 1. Standard LSC — `lsc_standard.m`

```
A = C + D = R' R ,       alpha = R \ (R' \ y)
```

* `O(N^2)` storage for `A`, `O(N^3)` work.
* The matrix is assembled block by block into a pre-allocated array
  (`build_covariance_matrix.m`); only the upper triangle is evaluated. A naive
  assembly would materialise a double-precision copy of the whole matrix first
  (20 GB instead of 10 GB at N = 5 × 10⁴).
* The default precision is **single**, 4 N² bytes. This is what makes
  N = 5 × 10⁴ fit into 32 GB; the coefficients are then accurate to ~1e-6
  relative, whereas a height anomaly in cm needs ~1e-2.
* `info.bytes = 4 N^2`.

## 2. Sequential LSC — `lsc_sequential.m`

Observations are sorted into spatially contiguous groups of `n_b`
(`order = 'lat'`: latitude strips), and processed recursively:

```
v_i     = y_i - A_{i,<i} alpha_{<i}
S_i     = A_ii - A_{i,<i} P_{<i} A_{<i,i}          (n_b x n_b, the only inverse)
K_i     = P_{<i} A_{<i,i} S_i^{-1}
alpha_{<i} <- alpha_{<i} - K_i v_i
alpha_i     = S_i^{-1} v_i
P_{<=i} = [ P_{<i} + K_i S_i K_i'   -K_i    ]
          [ -K_i'                   S_i^{-1}]
```

* The only matrix inverted is `n_b × n_b`, but the state `P` is `N × N`, so the
  **total storage is unchanged**: `info.pBytes = 8 N^2` (double) or `4 N^2`.
* Two products are evaluated row block by row block, so that the peak memory
  stays close to the size of `P` itself:
  * `P + K S K' = P + (P A) S^{-1} (P A)'`, accumulated in blocks of 2048 rows.
    The naive `K*S*K'` creates a `prev × prev` temporary (16.9 GB in double at
    N = 5 × 10⁴) which, on top of `P`, exceeds the memory of the machine. This
    identity is mentioned in Section 3.2 of the paper.
  * `PA = P(1:prev,1:prev) A'`, accumulated in blocks of 1024 rows. Writing it
    directly copies the whole leading block of `P` — as much memory as `P`
    again — which would double the peak memory of the recursion and make the
    paper's 18.6 GB figure unreachable in practice.
* Set `opt.single = true` (or `GSLSC_SEQ_SINGLE=1` for the script) to hold `P`
  in single precision: `4 N^2` bytes instead of `8 N^2`, at a relative accuracy
  of about 1e-6 in the coefficients.
* The solution is exact after one pass: no iteration and no stopping criterion.

## 3. Block Gauss-Seidel LSC — `lsc_block_gauss_seidel.m`

```
alpha_k^(m+1) = A_kk^{-1} ( y_k - sum_{j<k} A_kj alpha_j^(m+1)
                                 - sum_{j>k} A_kj alpha_j^(m) )
```

* Only the `n_b × n_b` diagonal blocks are factored; storage is
  `8 n_b^2 + O(N)` — 32 MB for `n_b = 2000`, 128 MB for `n_b = 4000` —
  **independent of N**.
* Every sweep evaluates all `N^2` covariance elements again. At N = 5 × 10⁴ one
  sweep is about 28 s, against 0.17 s for a matrix–vector product with the
  assembled matrix: the method trades time for storage.
* No state is carried across sweeps: the only unknown is `alpha`.

### Stopping rule

The coefficients converge much more slowly than the predicted quantity, so the
iteration stops on the **prediction** (Section 4.1 of the paper):

```
theta(alpha) = the height anomaly at the check points
stop when  max | theta(alpha^(m+1)) - theta(alpha^(m)) |  <  0.01 cm
```

`opt.thetaFun` supplies `theta`; `run_block_gs_lsc.m` builds it from
`signal_predict` and the 222 GNSS/levelling points. With this rule the Colorado
sub-region needs 33–40 sweeps at a residual of 8 × 10⁻³ – 2 × 10⁻². A residual
of 1 × 10⁻⁶ would need more than 500 sweeps for an improvement of less than
0.01 cm.

The residual reported in `info.relResid` is the exact
`||A alpha - y|| / ||y||`, computed with one extra matrix–vector product at the
end of the iteration.

## 4. Reference solver — `lsc_pcg.m`

Jacobi-preconditioned conjugate gradients on the same system. One iteration
costs one `lsc_matvec` (N² covariance elements); the residual decreases
monotonically. It is included for the comparison of Section 4, not as a
competitor of the three solvers: the number of iterations is problem dependent.

## 5. Complexity summary (N = 50 902)

| solver | matrix inverted | storage | work |
|---|---|---|---|
| standard | N × N | 4 N² = 9.65 GB | one Cholesky factorisation |
| sequential | n_b × n_b | 8 N² = 19.3 GB | one pass, N × (N − n_b) / 2 covariance blocks |
| block Gauss-Seidel | n_b × n_b | 8 n_b² = 32 MB / 128 MB | 33–40 full sweeps, N² elements each |

The storage figures quoted in the paper (9.31 GB, 18.63 GB, 32 MB, 128 MB) are
those of N = 5 × 10⁴; `storage_model.m` reproduces them.
