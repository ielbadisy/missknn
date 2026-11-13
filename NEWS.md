# missknn 1.2.2

* Multiple imputation (`m >= 2`) is now proper: each imputation draws its donors from an approximate Bayesian bootstrap (Rubin and Schenker 1986) of each variable's observed cases before the neighbor search and donor draw. Previously only the donor draw varied, so pooled intervals undercovered (83.5% in a 200-replicate simulation with a MAR covariate; 89.5% now, the rest due to kNN bias under that MAR mechanism, with Rubin's variance now close to the empirical variance). Multiple imputation results change; single imputation results do not. No API change.
* `numeric_estimator = "mean"` is now respected by the per-column tuning, which could previously still select local regression.
* The vignette and README pool with `mimar::poollm()` (the current name of the deprecated `mimar::pool_lm()`).

# missknn 1.2.1

* Internal call `basetable::as_basetable()` replaced by `basetable::asbasetable()` (the old name is deprecated in basetable 1.4.5); the minimum basetable version is now 1.4.5. No change to missknn's own API.

# missknn 1.2.0

## Breaking changes

* `missknn()` is now an imputer only and returns completed data directly: a
  `data.frame` when `m = 1`, and a list of `m` completed data frames with class
  `c("mimar_mi", "list")` when `m >= 2`.
* `complete()`, the `missknn` object class and its `print()`/`summary()`
  methods are removed. Analysis and pooling of multiply imputed data belong to
  dedicated packages such as `mimar` (`mimar::pool_lm()`).
* Completed data are plain data frames; `data.table` is no longer a dependency.
* Tabular handling (input coercion, column classes, type restoration) now
  uses `basetable` (>= 1.4.2). Only `NA` cells are imputed: blank strings stay
  a valid category.
* The `"tuning"` attribute reports the number of `missing` values per variable.
* Ordered factors are restored as ordered factors.

# missknn 1.1.2

* Base-R progress bar (`progress = TRUE`) with no extra dependency.
