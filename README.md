# rsomics-fligner

The Fligner-Killeen test for equality of variances — a value-exact, faster
drop-in for `scipy.stats.fligner`. It is the non-parametric, rank-based member
of the variance-homogeneity family (alongside Levene's and Bartlett's tests):
each group is centered, the pooled absolute deviations are jointly tie-averaged
ranked, the ranks are turned into normal scores, and the between-group score
dispersion is compared to a χ² with `k − 1` degrees of freedom. The
median-centered default is distribution-free and robust to non-normality.

## Install

```sh
cargo install rsomics-fligner
```

## Usage

Two or more single-column files, one group per file:

```sh
rsomics-fligner group1.tsv group2.tsv group3.tsv
```

Or one long-format file of `value<TAB>group` rows:

```sh
rsomics-fligner --long all.tsv
```

Output is a single tab-separated line `statistic<TAB>p`.

| Flag | Meaning |
|---|---|
| `--center {median,mean,trimmed}` | Centering statistic per group (default `median`). |
| `--proportiontocut <F>` | Fraction trimmed from each end when `--center trimmed` (default `0.05`). |
| `--long` | Treat the single input file as `value<TAB>group` rows. |
| `-t, --threads <N>` | Thread budget. |
| `--json` | Emit a single JSON result envelope. |
| `-q, --quiet` | Suppress progress on stderr. |

## Value-exactness

The statistic is computed in the same arithmetic order as
`scipy.stats.fligner` (`_morestats.py`), including NumPy's pairwise summation
for `np.mean` / `np.var(ddof=1)`. The two special functions are direct Cephes
ports rather than a third-party numerics crate:

- the rank→score transform `Φ⁻¹(rank/(2(N+1)) + 0.5)` uses a port of Cephes
  `ndtri` (= `scipy.special.ndtri` = `scipy.stats.norm.ppf`);
- the p-value `chi2.sf(statistic, k−1)` uses a port of Cephes `igamc` (=
  `scipy.special.chdtrc`).

Both are unit-tested against SciPy to ≤ 1e-12, and the committed
`tests/golden/` reproduce `scipy.stats.fligner` to ≤ 1e-13 (statistic) and
≤ 1e-12 (p-value) with no SciPy at test time.

## Origin

This crate is an independent Rust reimplementation of `scipy.stats.fligner`
based on:

- Fligner, M. A. and Killeen, T. J. (1976). *Distribution-free two-sample tests
  for scale.* Journal of the American Statistical Association 71(353), 210–213.
- Conover, W. J., Johnson, M. E. and Johnson, M. M. (1981). *A comparative study
  of tests for homogeneity of variances, with applications to the outer
  continental shelf bidding data.* Technometrics 23(4), 351–361
  (the F-K statistic, eq. 2.1, p. 355).
- The SciPy source `scipy/stats/_morestats.py::fligner` (scipy 1.17.1,
  BSD-3-Clause) for the exact transform and arithmetic order.
- The Cephes Math Library (Stephen L. Moshier) `ndtri.c` and `igam.c` for the
  inverse-normal and upper-incomplete-gamma special functions, the same code
  paths SciPy's `special` module uses.

License: MIT OR Apache-2.0.

Upstream credit: SciPy (<https://scipy.org>, BSD-3-Clause); Cephes
(<https://netlib.org/cephes/>, Moshier).
