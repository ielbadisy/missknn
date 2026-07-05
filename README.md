# missknn

`missknn` is an R package for fast masked k-nearest neighbor imputation of missing values in tabular data.

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
- multiple imputation can run multicore by default on Unix-like systems

It supports:

- single imputation
- multiple imputation
- numeric and mixed tabular data
- `data.table`-friendly workflows

## Benchmark

Two benchmark scripts write reports, plots, and CSV tables to `inst/benchmark/output/`:

```r
# accuracy/runtime vs. missForest, missMDA, VIM::kNN at n = 200, 30 simulations
Rscript inst/benchmark/benchmark_30.R

# runtime/accuracy vs. missMDA and missForest as n scales from 1,000 to 100,000
Rscript inst/benchmark/benchmark_bigdata.R
```

## Example

```r
library(missknn)

x <- data.frame(
  a = c(1, 2, NA, 4),
  b = c(2, NA, 3, 5),
  g = factor(c("x", "x", "y", NA))
)

imp <- missknn(x, k = 2, m = 1)
complete(imp)
```
