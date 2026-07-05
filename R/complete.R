#' Complete a missknn object
#'
#' @param object A `missknn` object.
#' @param action Either `"all"` to return all completed datasets or `"single"`
#'   to return one completed dataset. When `m = 1`, both options return the
#'   same dataset.
#' @param index Integer index of the completed dataset to return when
#'   `action = "single"` and `m > 1`.
#' @param ... Unused.
#'
#' @return A completed data.table or a list of completed data.tables.
#' @examples
#' x <- data.frame(a = c(1, NA, 3), b = c(2, 4, NA))
#' imp <- missknn(x, k = 1)
#' complete(imp)
#' @export
complete <- function(object, ...) {
  UseMethod("complete")
}

#' @rdname complete
#' @export
complete.missknn <- function(object, action = c("all", "single"), index = 1L, ...) {
  action <- match.arg(action)
  if (object$m == 1L) {
    return(object$completed)
  }
  if (action == "all") {
    return(object$completed)
  }
  if (index < 1L || index > length(object$completed)) {
    stop("`index` is out of range.", call. = FALSE)
  }
  object$completed[[index]]
}
