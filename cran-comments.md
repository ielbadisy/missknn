## Submission

This is a resubmission. missknn 1.0.0 is the currently published CRAN version;
this submission updates to 1.1.2.

## Changes since the CRAN release (1.0.0)

* Added a `progress` argument to `missknn()`, showing a live progress bar
  over the per-column masked-distance search on large datasets: one bar
  per tuning stage and per imputed column, updating in place with a
  percentage and an ETA. Implemented entirely in base R (no new package
  dependency).
* Requesting `progress = TRUE` runs the affected column's search on a
  single thread rather than in parallel (an explicit, documented
  visible-vs-fast tradeoff); `Rcpp::checkUserInterrupt()` provides Ctrl-C
  support. Automatically disabled when `m > 1` with `parallel_cores > 1`,
  since console output from `parallel::mclapply`'s forked worker processes
  cannot be shown live.

## Test environments

* local Ubuntu 24.04, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 1 note

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'
  (this flag comes from the local R toolchain configuration, not from the
  package's Makevars, and is not expected to appear on CRAN's build servers)

## Downstream dependencies

There are currently no downstream dependencies for this package.
