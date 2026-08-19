# missknn

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/ielbadisy/missknn)

`missknn` is an R package for fast masked k-nearest neighbor imputation of missing values in tabular data.

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("ielbadisy/missknn")
```

It is fast by design:

- distance computation, donor ranking, and aggregation are done in C++, with partial (not full)
  top-`k` selection and a batched per-column API that avoids per-cell R/C++ round trips
- each column's `k` and estimator (distance-weighted mean vs. a local weighted-regression fit)
  are chosen automatically from a holdout evaluation, rather than fixed globally
- columns with no locally exploitable signal are detected and filled directly from the global
  mean/mode, skipping neighbor search entirely
- a `donor_cap` bounds the neighbor-search candidate pool so cost stays roughly linear in `n`
  instead of quadratic, which matters at large `n`
- `data.table` is used for lightweight tabular orchestration
- multiple imputation (`m > 1`) can run multicore by default on Unix-like systems via `parallel::mclapply`
- within a single imputation, the per-column holdout `k`/estimator tuning and the deterministic
  (non-stochastic) neighbor search are parallelized across columns/receivers with `RcppParallel`,
  on all platforms including Windows

It supports:

- single imputation
- multiple imputation
- numeric and mixed tabular data
- `data.table`-friendly workflows

## About

`missknn` is developed and maintained by Imad El Badisy (elbadisyimad@gmail.com), released
under the MIT license.

## Latest Implementation

The current implementation parallelizes the per-column `k`/estimator tuning and the deterministic
single-imputation path with `RcppParallel`, including on Windows. The `m > 1` stochastic path still
uses process-level parallelism via `parallel::mclapply` on Unix-like systems.

## Parallel computing

`missknn` parallelizes at two independent levels:

- **Across imputations** (`m > 1`): each of the `m` completed datasets is generated in a
  separate process via `parallel::mclapply`, controlled by `parallel_cores` (defaults to
  `detectCores() - 1`). This only applies on Unix-like systems; it falls back to serial on
  Windows.
- **Within a single imputation**: the per-column holdout search that picks each column's `k`
  and estimator, and the deterministic (`m = 1`) neighbor search/aggregation, are split across
  threads with `RcppParallel::parallelFor` - across target columns for tuning, across receiver
  rows for imputation. This runs on every platform, including Windows, and benefits `m = 1` runs
  that get no benefit from `parallel_cores`. The stochastic sampling path used when `m > 1`
  stays single-threaded per process, since R's RNG isn't thread-safe; that case still parallelizes
  across the `m` processes instead.

Thread count for the `RcppParallel` layer follows the usual `RcppParallel::setThreadOptions()` /
`RCPP_PARALLEL_NUM_THREADS` conventions.

## Progress reporting

```r
imp <- missknn(big_data, k = 5, progress = TRUE)
```

With `progress = TRUE`, `missknn()` shows a live, self-contained (base R only, no extra
dependency) progress bar for each tuning stage (`Tuning numeric columns`,
`Tuning categorical columns`) and for each column that goes through the neighbor search
(`Imputing <column>`), updating in place with a percentage and an ETA, e.g.:

```
Imputing x2              [===========>------------------]  38%  ETA: 4s
```

Columns with no locally exploitable signal (filled directly from the global mean/mode, see the
`missknn()` API entry above) don't get a bar, since they skip the neighbor search entirely.

Because `RcppParallel` worker threads can't safely drive a progress display, requesting
`progress = TRUE` runs that column's search on a single thread instead of in parallel - an
explicit, visible-vs-fast tradeoff, not a bug. Leave `progress = FALSE` (the default) for the
fastest run. `progress` is also automatically disabled whenever `m > 1` with `parallel_cores > 1`,
since console output from `parallel::mclapply`'s forked worker processes can't be shown live.

## Benchmark

Two benchmark scripts write reports, plots, and CSV tables to `inst/benchmark/output/`:

```r
# accuracy/runtime vs. missForest, missMDA, VIM::kNN at n = 200, 30 simulations
Rscript inst/benchmark/benchmark_30.R

# runtime/accuracy vs. missMDA and missForest as n scales from 1,000 to 100,000
Rscript inst/benchmark/benchmark_bigdata.R
```

## API

### `missknn(data, k = 5L, m = 1L, scale = TRUE, add_indicator = FALSE, seed = NULL, weights = c("distance", "uniform"), numeric_estimator = c("regression", "mean"), ridge = 1e-4, donor_cap = 2000L, max_iter = 1L, parallel_cores = missknn_default_cores(), progress = FALSE)`

Fits masked KNN imputation and returns a `missknn` object.

| Argument | Description |
| --- | --- |
| `data` | A `data.frame` or matrix. |
| `k` | Number of nearest donors used per target column. |
| `m` | Number of completed datasets. `m = 1` is single (deterministic) imputation; `m > 1` is multiple (stochastic donor-sampling) imputation. |
| `scale` | Standardize numeric variables before computing distances. |
| `add_indicator` | Append `TRUE`/`FALSE` missingness-indicator columns to the completed output. |
| `seed` | Integer seed for reproducible multiple imputation. |
| `weights` | `"distance"` (inverse-distance weighted aggregation) or `"uniform"` (equal weights) over the `k` donors. |
| `numeric_estimator` | `"regression"` fits a distance-weighted local linear regression of the target on the receiver's observed numeric predictors over the `k` donors, falling back to the weighted mean when the local design is degenerate; `"mean"` always uses the weighted mean. |
| `ridge` | Ridge penalty added to the predictor covariance under `numeric_estimator = "regression"`. |
| `donor_cap` | Caps the donor pool searched per target column, subsampling once per column when exceeded, so search cost stays roughly linear in `n`. Does not affect the global mean/mode fallback, which always uses every donor. |
| `max_iter` | Number of single-pass refinement iterations over the whole dataset. |
| `parallel_cores` | Cores used across the `m` completed datasets when `m > 1` (via `parallel::mclapply`, Unix-like only). Defaults to `detectCores() - 1`. Independent of the `RcppParallel` threading used within a single imputation (see [Parallel computing](#parallel-computing)). |
| `progress` | Show a live `cli` progress bar per tuning stage and per imputed column (see [Progress reporting](#progress-reporting)). Default `FALSE`. |

Columns with no locally exploitable signal are detected during fitting and filled from the
global mean/mode directly, without a neighbor search.

Returns a `missknn` object: a list with `data` (input as a `data.table`), `completed` (one
completed data.table, or a list of `m` when `m > 1`), `imputations` (always a list of `m`),
plus the arguments used to fit (`k`, `m`, `scale`, `weights`, etc.) and internal `meta`.

### `complete(object, action = c("all", "single"), index = 1L, ...)`

S3 generic (`complete.missknn`) that extracts completed data from a `missknn` object.

- `action = "all"` (default): returns every completed dataset (a single data.table when `m = 1`,
  a list of `m` data.tables when `m > 1`).
- `action = "single"`: returns one completed data.table, selected by `index` (ignored when
  `m = 1`, since there is only one).

### `print(x)` / `summary(object)`

`print.missknn` reports `n`, `p`, `k`, `m`, `scale`, and `add_indicator`. `summary.missknn`
additionally reports per-variable and total missing counts (`missing_by_variable`,
`missing_total`) and has its own `print.summary.missknn` method.

## Example

```r
library(missknn)

x <- data.frame(
  a = c(1, 2, NA, 4),
  b = c(2, NA, 3, 5),
  g = factor(c("x", "x", "y", NA))
)

imp <- missknn(x, k = 2, m = 1)
imp
summary(imp)
complete(imp)

# multiple imputation
mi <- missknn(x, k = 2, m = 5, seed = 1)
complete(mi, action = "single", index = 2)
```

## Execution time vs. sample size

Runtime and NRMSE from n = 1,000 up to n = 100,000 (5 numeric columns, 20% MCAR per column).
`missForest` and `VIM::kNN` are only run up to n = 5,000 (per-tree cost and O(n^2) distance
search respectively make larger n impractical); `missknn` and `missRanger` run across the full
grid. Reproduce with:

```r
Rscript inst/benchmark/benchmark_bigdata.R
```

![Runtime plot](inst/benchmark/output/bigdata_runtime_plot.png)

![NRMSE plot](inst/benchmark/output/bigdata_mse_plot.png)

Table: Runtime and NRMSE by method and sample size

|      n|method     | runtime_sec|   nrmse|
|------:|:----------|-----------:|-------:|
|  1,000|missknn    |      0.0131| 1.00223|
|  1,000|missRanger |      0.1630| 1.15781|
|  1,000|missForest |      3.1537| 1.07461|
|  1,000|VIM::kNN   |      0.6088| 1.16132|
|  5,000|missknn    |      0.0321| 1.00020|
|  5,000|missRanger |      1.3632| 1.14478|
|  5,000|missForest |     75.1881| 1.03685|
|  5,000|VIM::kNN   |     10.9076| 1.15992|
| 20,000|missknn    |      0.0456| 0.99998|
| 20,000|missRanger |     10.3372| 1.14941|
|100,000|missknn    |      0.0938| 1.00000|
|100,000|missRanger |     85.6539| 1.14435|

`missForest`'s runtime already becomes impractical at n = 5,000 (over a minute) and would take
hours at n = 100,000; `VIM::kNN`'s distance search is similarly infeasible past a few thousand
rows. `missknn`'s capped-donor masked-KNN search stays close to linear in `n` and keeps running
in well under a second even at n = 100,000. `missRanger` also scales to n = 100,000 but at
substantially higher runtime.
