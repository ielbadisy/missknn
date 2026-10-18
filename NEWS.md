# missknn 1.2.0

## Breaking changes

* `missknn()` is now an imputer only and returns completed data directly: a
  `data.frame` when `m = 1`, and a list of `m` completed data frames with class
  `c("mimar_mi", "list")` when `m >= 2`.
* `complete()`, the `missknn` object class and its `print()`/`summary()`
  methods are removed. Analysis and pooling of multiply imputed data belong to
  dedicated packages such as `mimar` (`mimar::pool()`).
* Completed data are plain data frames; `data.table` is no longer a dependency.
* Tabular handling (input coercion, column classes, type restoration) now
  uses `basetable` (>= 1.4.2). Only `NA` cells are imputed: blank strings stay
  a valid category.
* The `"tuning"` attribute reports the number of `missing` values per variable.
* Ordered factors are restored as ordered factors.

# missknn 1.1.2

* Base-R progress bar (`progress = TRUE`) with no extra dependency.
