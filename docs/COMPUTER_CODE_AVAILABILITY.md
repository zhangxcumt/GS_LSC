# Computer Code Availability — text to paste into the manuscript

*Computers & Geosciences* requires a section at the end of the manuscript
describing the distribution of the code. Paste the section below into
`manuscript/20260923 GS-LSC (EN).docx`, after the "Data availability"
statement (before the declarations and the references; check the placement
required by the current CAGEO Word template).

**Two placeholders must be filled in before resubmission:**

1. `<GITHUB-ACCOUNT>` — the account that hosts the repository.
2. `<ZENODO-DOI>` — the DOI returned by Zenodo when the repository is archived
   (optional but strongly recommended; if you do not archive, delete the DOI
   clause and the DOI line in the README).

---

## Text for the manuscript

### Computer Code Availability

The code used in this study is open source and is distributed from a public
repository; it is **not** included as supplementary material.

| | |
|---|---|
| Name of the code | GS-LSC |
| Developers | Xun Zhang, Lubin Chang, Guobin Chang |
| Contact | Lubin Chang, 1410021052@nue.edu.cn |
| Year first available | 2026 |
| Hardware required | Any x86-64 workstation. The quick test and all examples run in under 8 GB of RAM; reproducing the full experiments of the paper (N = 50 902 observations) needs about 32 GB for the standard and sequential solvers, but only 128 MB for the block Gauss-Seidel solver proposed in the paper. |
| Program language | MATLAB |
| Program size | 35 MATLAB files, about 2 350 lines |
| Software required | MATLAB R2018b or newer; no toolbox |
| Availability | https://github.com/&lt;GITHUB-ACCOUNT&gt;/GS-LSC |
| License | MIT |
| DOI | &lt;ZENODO-DOI&gt; |

The repository contains the three solvers compared in the paper — standard
collocation with a Cholesky factorisation of the N × N normal-equation matrix,
sequential collocation in the sense of Moritz (1973) and Tscherning (1974), and
the block Gauss-Seidel iteration — together with the Tscherning–Rapp covariance
model, the Colorado test field of the 1-cm geoid experiment, the scripts that
reproduce Tables 3–4 and Figures 1, 4 and 5, and a README file with the
instructions to run them.

A self-contained quick test is provided in `examples/run_quick_test.m`. It needs
no input data beyond that distributed with the repository and runs in about one
minute: it first checks the three solvers against a known analytical solution,
and then compares their height anomaly predictions at the 222 GNSS/levelling
check points of the Colorado sub-region. Both checks print their tolerances and
stop with an error if a check fails. The README file gives the commands to run
the quick test and to reproduce each table and figure of the paper.

---

## Checklist against the journal's requirements

| Requirement (from the editorial letter) | Where it is met |
|---|---|
| Section "Computer Code Availability" at the end of the manuscript | this document |
| Public repository, anonymous download | GitHub public repository, no login required |
| Not a single compressed file; individual source files | repository holds 35 individual `.m` files |
| README with the purpose of the code and how to use it | `README.md` |
| At least one quick test or example | `examples/example_known_solution.m`, `examples/example_colorado_subset.m` |
| Instructions in the README on how to run the test | `README.md`, Section 3 |
| Open-source license (optional but recommended) | `LICENSE` (MIT) |
| No proprietary software needed for the presented work | MATLAB only; no toolbox; all computations are in the repository |
| Code not uploaded as supplementary material | remove any code file from the submission files in Editorial Manager |
