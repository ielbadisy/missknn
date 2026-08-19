missknn_default_cores <- function() {
  # CRAN policy caps examples/tests/vignettes at 2 cores; check-time runs set
  # `_R_CHECK_LIMIT_CORES_`, in which case fall back to that limit instead of
  # `detectCores() - 1`, which can otherwise try to spawn far more workers
  # than the check environment allows.
  limited <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  if (!is.na(limited) && nzchar(limited)) {
    return(2L)
  }
  max(1L, parallel::detectCores() - 1L)
}

# A cli bar is opened per column/tuning-batch (labelled by name) and closed
# right after that call returns; `.auto_close = FALSE` because relying on
# cli's default frame-based auto-close would only clear the bar when the
# *enclosing* function returns, not after each column, leaving finished
# bars stacked on screen instead of replaced in place.
missknn_progress_start <- function(label, total, enabled) {
  if (!isTRUE(enabled) || total <= 0) {
    return(NULL)
  }
  cli::cli_progress_bar(name = label, total = total, clear = FALSE, .auto_close = FALSE)
}

missknn_progress_cb <- function(id) {
  if (is.null(id)) {
    return(NULL)
  }
  function(done, total) cli::cli_progress_update(id = id, set = done)
}

missknn_progress_done <- function(id) {
  if (!is.null(id)) {
    cli::cli_progress_done(id = id)
  }
  invisible(NULL)
}

missknn_validate_input <- function(data) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.", call. = FALSE)
  }

  dt <- data.table::as.data.table(data)
  data.table::copy(dt)
}

missknn_detect_types <- function(dt) {
  classes <- vapply(dt, function(x) {
    if (is.numeric(x)) {
      "numeric"
    } else if (is.logical(x) || is.factor(x) || is.character(x)) {
      "categorical"
    } else {
      stop("Unsupported column type: ", paste(class(x), collapse = "/"), call. = FALSE)
    }
  }, character(1))

  list(
    classes = classes,
    numeric_idx = which(classes == "numeric"),
    categorical_idx = which(classes == "categorical"),
    names = names(dt)
  )
}

missknn_prepare_numeric <- function(dt, numeric_idx, scale = TRUE) {
  n <- nrow(dt)
  p <- length(numeric_idx)
  if (!p) {
    return(list(
      scaled_matrix = matrix(numeric(0), nrow = n, ncol = 0),
      raw_matrix = matrix(numeric(0), nrow = n, ncol = 0),
      center = numeric(0),
      scale = numeric(0)
    ))
  }

  raw_mat <- matrix(NA_real_, nrow = n, ncol = p)
  scaled_mat <- matrix(NA_real_, nrow = n, ncol = p)
  center <- numeric(p)
  scale_vec <- rep(1, p)

  for (k in seq_along(numeric_idx)) {
    x <- as.numeric(dt[[numeric_idx[k]]])
    raw_mat[, k] <- x
    obs <- x[!is.na(x)]
    center[k] <- if (length(obs)) mean(obs) else 0
    s <- if (length(obs) > 1) stats::sd(obs) else 0
    scale_vec[k] <- if (isTRUE(scale) && is.finite(s) && s > 0) s else 1
    if (isTRUE(scale)) {
      scaled_mat[, k] <- (raw_mat[, k] - center[k]) / scale_vec[k]
    } else {
      scaled_mat[, k] <- raw_mat[, k]
    }
  }

  list(scaled_matrix = scaled_mat, raw_matrix = raw_mat, center = center, scale = scale_vec)
}

missknn_prepare_categorical <- function(dt, categorical_idx) {
  n <- nrow(dt)
  p <- length(categorical_idx)

  if (!p) {
    return(list(matrix = matrix(integer(0), nrow = n, ncol = 0), levels = list(), original_types = character(0)))
  }

  mat <- matrix(NA_integer_, nrow = n, ncol = p)
  levels_list <- vector("list", p)
  original_types <- character(p)

  for (k in seq_along(categorical_idx)) {
    x <- dt[[categorical_idx[k]]]
    original_types[k] <- paste(class(x), collapse = "/")
    if (is.logical(x)) {
      lev <- c("FALSE", "TRUE")
      mat[, k] <- match(as.character(x), lev)
      levels_list[[k]] <- lev
    } else if (is.factor(x)) {
      lev <- levels(x)
      mat[, k] <- as.integer(x)
      levels_list[[k]] <- lev
    } else {
      x_chr <- as.character(x)
      lev <- sort(unique(stats::na.omit(x_chr)))
      mat[, k] <- match(x_chr, lev)
      levels_list[[k]] <- lev
    }
  }

  list(matrix = mat, levels = levels_list, original_types = original_types)
}

missknn_mode <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

missknn_weighted_mode <- function(values, weights) {
  keep <- !is.na(values) & !is.na(weights)
  values <- values[keep]
  weights <- weights[keep]
  if (!length(values)) {
    return(NA_character_)
  }
  agg <- tapply(weights, values, sum)
  agg <- sort(agg, decreasing = TRUE)
  names(agg)[1]
}

missknn_global_fill <- function(target, meta) {
  if (meta$classes[target] == "numeric") {
    obs <- meta$original[[target]]
    obs <- obs[!is.na(obs)]
    if (!length(obs)) return(NA_real_)
    mean(obs)
  } else {
    vals <- meta$original[[target]]
    vals <- vals[!is.na(vals)]
    if (!length(vals)) return(NA_character_)
    missknn_mode(as.character(vals))
  }
}

missknn_decode_categorical <- function(code, meta, target) {
  if (is.na(code)) {
    return(NA_character_)
  }
  lev <- meta$categorical$levels[[meta$categorical_pos[target]]]
  as.character(lev[code])
}

missknn_distance <- function(i, donors, target, working, meta) {
  if (!length(donors)) {
    return(numeric(0))
  }

  sum_w <- numeric(length(donors))
  sum_d <- numeric(length(donors))

  for (col in seq_along(meta$classes)) {
    if (col == target) {
      next
    }
    wj <- meta$weights[col]
    if (wj <= 0) {
      next
    }

    if (meta$classes[col] == "numeric") {
      r <- working$numeric[i, meta$numeric_pos[col]]
      d <- working$numeric[donors, meta$numeric_pos[col]]
      mask <- !is.na(r) & !is.na(d)
      if (any(mask)) {
        diff2 <- (d[mask] - r)^2
        sum_d[mask] <- sum_d[mask] + wj * diff2
        sum_w[mask] <- sum_w[mask] + wj
      }
    } else {
      r <- working$categorical[i, meta$categorical_pos[col]]
      d <- working$categorical[donors, meta$categorical_pos[col]]
      mask <- !is.na(r) & !is.na(d)
      if (any(mask)) {
        diff <- as.numeric(d[mask] != r)
        sum_d[mask] <- sum_d[mask] + wj * diff
        sum_w[mask] <- sum_w[mask] + wj
      }
    }
  }

  out <- rep(Inf, length(donors))
  ok <- sum_w > 0
  out[ok] <- sqrt(sum_d[ok] / sum_w[ok])
  out
}

missknn_pick_donors <- function(distances, donors, k) {
  finite <- is.finite(distances)
  if (!any(finite)) {
    return(integer(0))
  }
  ord <- order(distances[finite], donors[finite])
  donor_sel <- donors[finite][ord]
  donor_sel[seq_len(min(k, length(donor_sel)))]
}

missknn_aggregate_single <- function(target, donor_rows, distances, working, meta) {
  if (!length(donor_rows)) {
    return(missknn_global_fill(target, meta))
  }

  donor_values <- working$original[donor_rows, target]
  d <- distances[match(donor_rows, meta$current_donor_index)]
  if (meta$aggregation == "uniform") {
    w <- rep(1, length(d))
  } else {
    w <- 1 / (d + meta$epsilon)
  }

  if (meta$classes[target] == "numeric") {
    if (all(!is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
      return(mean(as.numeric(donor_values), na.rm = TRUE))
    }
    stats::weighted.mean(as.numeric(donor_values), w = w, na.rm = TRUE)
  } else {
    keep <- !is.na(donor_values)
    donor_values <- donor_values[keep]
    w <- w[keep]
    if (!length(donor_values)) {
      return(missknn_global_fill(target, meta))
    }
    values_chr <- as.character(donor_values)
    agg <- tapply(w, values_chr, sum)
    agg <- sort(agg, decreasing = TRUE)
    names(agg)[1]
  }
}

missknn_sample_donor <- function(target, donor_rows, distances, working, meta) {
  if (!length(donor_rows)) {
    return(missknn_global_fill(target, meta))
  }

  d <- distances[match(donor_rows, meta$current_donor_index)]
  if (meta$aggregation == "uniform") {
    w <- rep(1, length(d))
  } else {
    w <- 1 / (d + meta$epsilon)
  }
  w[!is.finite(w)] <- 0
  if (sum(w) <= 0) {
    probs <- rep(1 / length(donor_rows), length(donor_rows))
  } else {
    probs <- w / sum(w)
  }

  pick <- sample(donor_rows, size = 1L, prob = probs)
  working$original[pick, target]
}

missknn_single_pass <- function(dt, meta, stochastic = FALSE) {
  # R's copy-on-modify semantics already duplicate these matrices the first
  # time `working$*[...] <-` writes into them, so no explicit copy is needed
  # up front; `meta`'s originals are never mutated in place.
  working <- list(
    original = as.data.frame(dt, stringsAsFactors = FALSE),
    numeric_scaled = meta$numeric$scaled_matrix,
    numeric_raw = meta$numeric$raw_matrix,
    categorical = meta$categorical$matrix
  )

  if (length(meta$categorical_idx)) {
    for (pos in seq_along(meta$categorical_idx)) {
      j <- meta$categorical_idx[pos]
      working$original[[j]] <- as.character(working$original[[j]])
    }
  }

  for (target in seq_along(meta$classes)) {
    receivers <- which(is.na(working$original[[target]]))
    if (!length(receivers)) {
      next
    }
    donors <- which(!is.na(working$original[[target]]))
    if (!length(donors)) {
      fill <- missknn_global_fill(target, meta)
      for (i in receivers) {
        working$original[[target]][i] <- fill
      }
      next
    }

    if (isTRUE(meta$is_global_col[target])) {
      # Validated by holdout tuning to carry no locally exploitable signal:
      # skip the O(receivers x donors) masked-distance search and fill
      # directly from the global donor mean/mode (or a bootstrap draw of it
      # for multiple imputation), which is both faster and at least as
      # accurate as any neighborhood-weighted estimate for this column.
      if (meta$classes[target] == "numeric") {
        donor_values <- working$numeric_raw[donors, meta$numeric_pos[target]]
        donor_values <- donor_values[!is.na(donor_values)]
        fills <- if (stochastic) {
          sample(donor_values, size = length(receivers), replace = TRUE)
        } else {
          rep(mean(donor_values), length(receivers))
        }
        working$original[[target]][receivers] <- fills
        working$numeric_raw[receivers, meta$numeric_pos[target]] <- fills
        if (isTRUE(meta$scale)) {
          working$numeric_scaled[receivers, meta$numeric_pos[target]] <-
            (fills - meta$numeric$center[meta$numeric_pos[target]]) / meta$numeric$scale[meta$numeric_pos[target]]
        } else {
          working$numeric_scaled[receivers, meta$numeric_pos[target]] <- fills
        }
      } else {
        donor_codes <- working$categorical[donors, meta$categorical_pos[target]]
        donor_codes <- donor_codes[!is.na(donor_codes)]
        codes <- if (stochastic) {
          sample(donor_codes, size = length(receivers), replace = TRUE)
        } else {
          tab <- table(donor_codes)
          rep(as.integer(names(tab)[which.max(tab)]), length(receivers))
        }
        fills <- vapply(codes, missknn_decode_categorical, character(1), meta = meta, target = target)
        working$original[[target]][receivers] <- fills
        working$categorical[receivers, meta$categorical_pos[target]] <- codes
      }
      next
    }

    # Bound the neighbor-search candidate pool so distance search stays
    # roughly linear in n instead of quadratic; the global mean/mode path
    # above already used the uncapped `donors` for its O(1)-per-receiver fill.
    donors_search <- if (length(donors) > meta$donor_cap) {
      sample(donors, meta$donor_cap)
    } else {
      donors
    }

    if (meta$classes[target] == "numeric") {
      col_k <- meta$k_per_col[target]
      col_estimator <- meta$estimator_per_col[target]
      # For "regression", k_per_col already stores the tuned regression
      # bandwidth directly (picked from the same adaptive grid used for the
      # weighted-mean k, spanning small local neighborhoods up to the entire
      # donor pool), so it is used as-is rather than rescaled again here.
      col_reg_neighbors <- min(col_k, length(donors_search))
      pb_id <- missknn_progress_start(
        paste0("Imputing ", names(meta$classes)[target]),
        length(receivers), isTRUE(meta$progress)
      )
      fills <- cpp_missknn_impute_numeric_column(
        working$numeric_scaled,
        working$numeric_raw,
        working$categorical,
        meta$type_codes,
        meta$numeric_pos,
        meta$categorical_pos,
        donors_search,
        receivers,
        target,
        meta$weights,
        col_k,
        meta$aggregation,
        col_estimator,
        stochastic,
        meta$epsilon,
        meta$ridge,
        col_reg_neighbors,
        isTRUE(meta$progress),
        missknn_progress_cb(pb_id)
      )
      missknn_progress_done(pb_id)
      working$original[[target]][receivers] <- fills
      working$numeric_raw[receivers, meta$numeric_pos[target]] <- fills
      if (isTRUE(meta$scale)) {
        working$numeric_scaled[receivers, meta$numeric_pos[target]] <-
          (fills - meta$numeric$center[meta$numeric_pos[target]]) / meta$numeric$scale[meta$numeric_pos[target]]
      } else {
        working$numeric_scaled[receivers, meta$numeric_pos[target]] <- fills
      }
    } else {
      pb_id <- missknn_progress_start(
        paste0("Imputing ", names(meta$classes)[target]),
        length(receivers), isTRUE(meta$progress)
      )
      codes <- cpp_missknn_impute_categorical_column(
        working$numeric_scaled,
        working$categorical,
        meta$type_codes,
        meta$numeric_pos,
        meta$categorical_pos,
        donors_search,
        receivers,
        target,
        meta$weights,
        meta$k_per_col[target],
        meta$aggregation,
        stochastic,
        meta$epsilon,
        isTRUE(meta$progress),
        missknn_progress_cb(pb_id)
      )
      missknn_progress_done(pb_id)
      fills <- vapply(codes, missknn_decode_categorical, character(1), meta = meta, target = target)
      working$original[[target]][receivers] <- fills
      working$categorical[receivers, meta$categorical_pos[target]] <- codes
    }
  }

  working$original
}

missknn_restore_types <- function(dt, meta) {
  out <- as.data.frame(dt, stringsAsFactors = FALSE)
  for (j in seq_along(meta$classes)) {
    if (meta$classes[j] == "numeric") {
      out[[j]] <- as.numeric(out[[j]])
    } else {
      lev <- meta$categorical$levels[[meta$categorical_pos[j]]]
      vals <- as.character(out[[j]])
      if (grepl("factor", meta$original_classes[j], fixed = TRUE)) {
        out[[j]] <- factor(vals, levels = lev)
      } else if (grepl("logical", meta$original_classes[j], fixed = TRUE)) {
        out[[j]] <- as.logical(vals)
      } else {
        out[[j]] <- vals
      }
    }
  }
  data.table::as.data.table(out)
}

missknn_candidate_k <- function(n_train, k_default) {
  # A single geometric bandwidth grid from small/local to the entire training
  # pool, shared by the numeric weighted-mean search, the local-regression
  # bandwidth search, and the categorical vote search. Spanning the full
  # range in one grid -- rather than a small local grid plus a separate
  # hardcoded "wide" special case -- means "how local should this column's
  # estimate be" is one adaptive choice along a continuum, not a discrete
  # local-vs-global switch bolted on afterwards. The O(1) global mean/mode
  # fast path (Section 2.4) is still evaluated analytically as a baseline
  # rather than via this grid, since it needs no distance search at all.
  base <- sort(unique(c(k_default, 5L, 15L, 40L)))
  base <- base[base <= n_train]
  if (!length(base)) {
    base <- max(1L, n_train)
  }
  as.integer(base)
}

# Data-driven k (and, for numeric columns, estimator) selection per column:
# holds out a subset of each column's observed values, imputes them from the
# remaining donors under each candidate (k, estimator), and keeps whichever
# minimizes holdout error. This lets a column shrink towards a large k (close
# to the global mean/mode) when it carries no local signal, or towards a small
# k when nearby donors are genuinely informative, instead of using one fixed k
# for every column.
missknn_tune_all <- function(dt, meta) {
  p <- length(meta$classes)
  k_per_col <- rep(meta$k, p)
  estimator_per_col <- rep(meta$numeric_estimator, p)
  is_global_col <- rep(FALSE, p)

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(20260705L)

  num_targets <- integer(0)
  num_train <- list()
  num_hold <- list()
  num_kgrid <- list()
  cat_targets <- integer(0)
  cat_train <- list()
  cat_hold <- list()
  cat_kgrid <- list()

  for (target in seq_len(p)) {
    donors <- which(!is.na(dt[[target]]))
    n_don <- length(donors)
    if (n_don < 20L) {
      next
    }

    n_hold <- max(5L, min(100L, floor(0.25 * n_don)))
    hold_idx <- sample(donors, n_hold)
    train_idx <- setdiff(donors, hold_idx)
    if (length(train_idx) < 5L) {
      next
    }
    if (length(train_idx) > meta$donor_cap) {
      train_idx <- sample(train_idx, meta$donor_cap)
    }

    candidates_k <- missknn_candidate_k(length(train_idx), meta$k)

    if (meta$classes[target] == "numeric") {
      num_targets <- c(num_targets, target)
      num_train[[length(num_train) + 1L]] <- train_idx
      num_hold[[length(num_hold) + 1L]] <- hold_idx
      num_kgrid[[length(num_kgrid) + 1L]] <- candidates_k
    } else {
      cat_targets <- c(cat_targets, target)
      cat_train[[length(cat_train) + 1L]] <- train_idx
      cat_hold[[length(cat_hold) + 1L]] <- hold_idx
      cat_kgrid[[length(cat_kgrid) + 1L]] <- candidates_k
    }
  }

  if (length(num_targets)) {
    pb_id <- missknn_progress_start("Tuning numeric columns", length(num_targets), isTRUE(meta$progress))
    res <- cpp_missknn_tune_numeric(
      meta$numeric$scaled_matrix, meta$numeric$raw_matrix, meta$categorical$matrix,
      meta$type_codes, meta$numeric_pos, meta$categorical_pos, meta$weights,
      meta$epsilon, meta$ridge, num_targets, num_train, num_hold, num_kgrid,
      isTRUE(meta$progress), missknn_progress_cb(pb_id)
    )
    missknn_progress_done(pb_id)
    k_per_col[num_targets] <- res$k
    estimator_per_col[num_targets] <- res$estimator
    is_global_col[num_targets] <- res$is_global
  }

  if (length(cat_targets)) {
    pb_id <- missknn_progress_start("Tuning categorical columns", length(cat_targets), isTRUE(meta$progress))
    res_cat <- cpp_missknn_tune_categorical(
      meta$numeric$scaled_matrix, meta$categorical$matrix,
      meta$type_codes, meta$numeric_pos, meta$categorical_pos, meta$weights,
      meta$epsilon, cat_targets, cat_train, cat_hold, cat_kgrid,
      isTRUE(meta$progress), missknn_progress_cb(pb_id)
    )
    missknn_progress_done(pb_id)
    k_per_col[cat_targets] <- res_cat$k
    is_global_col[cat_targets] <- res_cat$is_global
  }

  list(k_per_col = as.integer(k_per_col), estimator_per_col = estimator_per_col, is_global_col = is_global_col)
}

missknn_make_meta <- function(dt, k, scale, weights, add_indicator, seed,
                               numeric_estimator = "regression", ridge = 1e-4,
                               donor_cap = 2000L, progress = FALSE) {
  types <- missknn_detect_types(dt)
  reg_neighbors <- max(4L * as.integer(k), 30L)
  numeric_prep <- missknn_prepare_numeric(dt, types$numeric_idx, scale = scale)
  categorical_prep <- missknn_prepare_categorical(dt, types$categorical_idx)

  numeric_pos <- rep(0L, length(types$classes))
  categorical_pos <- rep(0L, length(types$classes))
  if (length(types$numeric_idx)) {
    numeric_pos[types$numeric_idx] <- seq_along(types$numeric_idx)
  }
  if (length(types$categorical_idx)) {
    categorical_pos[types$categorical_idx] <- seq_along(types$categorical_idx)
  }

  meta <- list(
    classes = types$classes,
    type_codes = ifelse(types$classes == "numeric", 1L, 2L),
    original_classes = vapply(dt, function(x) paste(class(x), collapse = "/"), character(1)),
    names = types$names,
    numeric_idx = types$numeric_idx,
    categorical_idx = types$categorical_idx,
    numeric_pos = numeric_pos,
    categorical_pos = categorical_pos,
    numeric = numeric_prep,
    categorical = categorical_prep,
    weights = rep(1, length(types$classes)),
    k = as.integer(k),
    scale = scale,
    add_indicator = add_indicator,
    aggregation = weights,
    numeric_estimator = numeric_estimator,
    ridge = ridge,
    donor_cap = as.integer(donor_cap),
    reg_neighbors = as.integer(reg_neighbors),
    epsilon = 1e-08,
    seed = seed,
    original = dt,
    current_donor_index = integer(0),
    progress = isTRUE(progress)
  )

  tuned <- missknn_tune_all(dt, meta)
  meta$k_per_col <- tuned$k_per_col
  meta$estimator_per_col <- tuned$estimator_per_col
  meta$is_global_col <- tuned$is_global_col
  meta
}

missknn_finalize <- function(completed, original, meta, add_indicator) {
  out <- missknn_restore_types(completed, meta)
  if (isTRUE(add_indicator)) {
    miss_cols <- names(original)[vapply(original, function(x) any(is.na(x)), logical(1))]
    for (nm in miss_cols) {
      out[[paste0(nm, "_missing")]] <- as.integer(is.na(original[[nm]]))
    }
  }
  out
}
