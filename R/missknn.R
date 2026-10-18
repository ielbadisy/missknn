#' Fast masked KNN imputation
#'
#' `missknn()` performs single imputation or multiple imputation for tabular
#' data using masked k-nearest neighbor distances computed only on jointly
#' observed variables.
#'
#' @param data A data.frame or matrix.
#' @param k Number of nearest donors to use.
#' @param m Number of completed datasets to generate. Use `m = 1` for single
#'   imputation.
#' @param scale Logical; whether to standardize numeric variables before
#'   computing distances.
#' @param add_indicator Logical; if `TRUE`, missingness indicator columns are
#'   appended to the completed output.
#' @param seed Optional integer seed for reproducible multiple imputation.
#' @param weights Character string. Currently `"distance"` uses inverse-distance
#'   aggregation and `"uniform"` uses equal weights.
#' @param numeric_estimator Character string. `"regression"` (default) fits a
#'   distance-weighted local linear regression of the target on the receiver's
#'   observed numeric predictors over the `k` nearest donors, falling back to
#'   the weighted mean when the local design is degenerate. `"mean"` always
#'   uses the distance-weighted mean.
#' @param ridge Ridge penalty added to the predictor covariance when
#'   `numeric_estimator = "regression"`.
#' @param donor_cap Maximum donor pool size searched per target column. When a
#'   column has more donors than this, a random subsample of `donor_cap`
#'   donors is drawn once per column and used as the neighbor-search
#'   candidate set, keeping distance-search cost roughly linear in `n` instead
#'   of quadratic. Does not affect the global mean/mode fast path, which
#'   always uses every donor.
#' @param max_iter Number of single-pass refinement iterations.
#' @param parallel_cores Number of cores used when `m > 1`. Defaults to the
#'   available cores minus one, capped at 2 under `R CMD check`-style
#'   environments that set `_R_CHECK_LIMIT_CORES_`.
#' @param progress Logical; if `TRUE`, shows a console progress bar over the
#'   receiver rows (or, during tuning, over target columns) of each column's
#'   masked-distance search, which is where compute time concentrates on
#'   large datasets. Automatically disabled when `m > 1` and imputations run
#'   in parallel forked processes (`parallel_cores > 1`), since console
#'   output from forked workers cannot be displayed live. Default `FALSE`.
#'
#' @details
#' Let \eqn{X = (x_{ij})} be an \eqn{n \times p} data matrix and let \eqn{R_{ij}} be the
#' missingness indicator:
#'
#' \deqn{R_{ij} = 1 \; \text{if } x_{ij} \text{ is observed}, \qquad
#'       R_{ij} = 0 \; \text{if } x_{ij} \text{ is missing}.}
#'
#' For a receiver row `i`, a donor row `l`, and a target variable `t`,
#' the masked distance is computed on the shared observed set
#'
#' \deqn{S_{ilt} = O_i \cap O_l \setminus \{t\}.}
#'
#' The mixed-data masked distance is
#'
#' \deqn{
#' d^2_t(i,l) =
#' \frac{\sum_{j \in S_{ilt}} w_j \, delta_j(i,l)}
#' {\sum_{j \in S_{ilt}} w_j},
#' }
#'
#' where `delta_j(i,l)` is squared difference for numeric variables and an
#' indicator of inequality for categorical variables.
#'
#' Single imputation uses a distance-weighted estimator from the `k` nearest
#' donors. Multiple imputation samples a donor from the same `k` nearest
#' donors using normalized inverse-distance probabilities.
#'
#' `missknn()` is an imputer only: it returns completed data. Analysis and
#' pooling of multiply imputed data (fitting a model on each completed data
#' set and combining the fits with Rubin's rules) are left to dedicated
#' packages such as `mimar` (for example `mimar::pool_lm()` on a list of
#' `lm()` fits). Tabular handling (input coercion, column classes, type
#' restoration) uses the `basetable` package. Only `NA` cells are imputed.
#'
#' @return With `m = 1`, a completed `data.frame` with the same columns and
#'   column types as `data` (plus `<name>_missing` indicator columns when
#'   `add_indicator = TRUE`). With `m >= 2`, a list of `m` such data frames
#'   with class `c("mimar_mi", "list")`. Either carries a `"tuning"`
#'   attribute: a data frame with one row per variable giving its number of
#'   `missing` values, the tuned number
#'   of donors `k`, the `estimator` (`"regression"`, `"mean"`, or `"vote"` for
#'   categorical variables) and whether the variable was filled from its
#'   global mean or mode (`global`, no neighbor signal found).
#' @examples
#' set.seed(1)
#' x <- data.frame(
#'   a = c(1, 2, NA, 4, 5),
#'   b = c(2, NA, 3, 4, 6),
#'   g = factor(c("x", "x", "y", NA, "y"))
#' )
#' missknn(x, k = 2)
#'
#' # multiple imputation: fit on each completed data set
#' imps <- missknn(x, k = 2, m = 3, seed = 1, parallel_cores = 1L)
#' fits <- lapply(imps, function(d) lm(a ~ b, data = d))
#' sapply(fits, coef)
#' @export
missknn <- function(data, k = 5L, m = 1L, scale = TRUE, add_indicator = FALSE,
                    seed = NULL, weights = c("distance", "uniform"),
                    numeric_estimator = c("regression", "mean"), ridge = 1e-4,
                    donor_cap = 2000L,
                    max_iter = 1L,
                    parallel_cores = missknn_default_cores(),
                    progress = FALSE) {
  weights <- match.arg(weights)
  numeric_estimator <- match.arg(numeric_estimator)
  dt <- missknn_validate_input(data)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (k < 1) stop("`k` must be at least 1.", call. = FALSE)
  if (m < 1) stop("`m` must be at least 1.", call. = FALSE)
  if (max_iter < 1) stop("`max_iter` must be at least 1.", call. = FALSE)
  if (donor_cap < 1) stop("`donor_cap` must be at least 1.", call. = FALSE)
  if (m > 1L && parallel_cores < 1) stop("`parallel_cores` must be at least 1.", call. = FALSE)

  runs_forked <- m > 1L && parallel_cores > 1L && .Platform$OS.type == "unix"
  use_progress <- isTRUE(progress) && !runs_forked

  meta <- missknn_make_meta(dt, k = k, scale = scale, weights = weights, add_indicator = add_indicator,
                             seed = seed, numeric_estimator = numeric_estimator, ridge = ridge,
                             donor_cap = donor_cap, progress = use_progress)

  run_one <- function(r) {
    # dt is a basetable (a data.frame), so copy-on-modify keeps it intact across runs.
    working <- dt
    for (iter in seq_len(max_iter)) {
      working <- missknn_single_pass(working, meta, stochastic = m > 1L)
    }
    missknn_finalize(working, dt, meta, add_indicator)
  }

  if (runs_forked) {
    completed_sets <- parallel::mclapply(seq_len(m), run_one, mc.cores = parallel_cores)
  } else {
    completed_sets <- lapply(seq_len(m), run_one)
  }

  tuning <- data.frame(
    variable = names(dt),
    missing = as.integer(meta$missing),
    k = as.integer(meta$k_per_col),
    estimator = ifelse(meta$classes == "numeric", meta$estimator_per_col, "vote"),
    global = as.logical(meta$is_global_col),
    stringsAsFactors = FALSE
  )
  out <- if (m == 1L) completed_sets[[1L]] else structure(completed_sets, class = c("mimar_mi", "list"))
  attr(out, "tuning") <- tuning
  out
}
