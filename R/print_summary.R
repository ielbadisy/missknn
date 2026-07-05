#' @rdname missknn
#' @export
print.missknn <- function(x, ...) {
  cat("<missknn>\n", sep = "")
  cat("  observations: ", nrow(x$data), "\n", sep = "")
  cat("  variables: ", ncol(x$data), "\n", sep = "")
  cat("  k: ", x$k, "\n", sep = "")
  cat("  completed datasets: ", x$m, "\n", sep = "")
  cat("  scale: ", x$scale, "\n", sep = "")
  cat("  add_indicator: ", x$add_indicator, "\n", sep = "")
  invisible(x)
}

#' @rdname missknn
#' @export
summary.missknn <- function(object, ...) {
  missings <- vapply(object$data, function(x) sum(is.na(x)), integer(1))
  out <- list(
    n = nrow(object$data),
    p = ncol(object$data),
    k = object$k,
    m = object$m,
    scale = object$scale,
    add_indicator = object$add_indicator,
    missing_by_variable = missings,
    missing_total = sum(missings)
  )
  class(out) <- "summary.missknn"
  out
}

#' @rdname missknn
#' @export
print.summary.missknn <- function(x, ...) {
  cat("<summary.missknn>\n", sep = "")
  cat("  n: ", x$n, "\n", sep = "")
  cat("  p: ", x$p, "\n", sep = "")
  cat("  k: ", x$k, "\n", sep = "")
  cat("  m: ", x$m, "\n", sep = "")
  cat("  scale: ", x$scale, "\n", sep = "")
  cat("  add_indicator: ", x$add_indicator, "\n", sep = "")
  cat("  total missing: ", x$missing_total, "\n", sep = "")
  invisible(x)
}
