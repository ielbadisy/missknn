mixed <- function() {
  set.seed(3)
  n <- 60
  x <- data.frame(
    a = rnorm(n),
    b = rnorm(n),
    g = factor(sample(c("u", "v", "w"), n, TRUE)),
    l = sample(c(TRUE, FALSE), n, TRUE),
    s = sample(c("p", "q"), n, TRUE),
    stringsAsFactors = FALSE
  )
  for (j in names(x)) x[[j]][sample(n, 8)] <- NA
  x
}

test_that("m = 1 returns a completed data.frame with the input types", {
  x <- mixed()
  out <- missknn(x, k = 3)
  expect_s3_class(out, "data.frame")
  expect_identical(class(out), "data.frame")
  expect_identical(names(out), names(x))
  expect_identical(dim(out), dim(x))
  expect_false(anyNA(out))
  expect_identical(lapply(out, class), lapply(x, class))
  expect_identical(levels(out$g), levels(x$g))
  obs <- !is.na(x$a)
  expect_equal(out$a[obs], x$a[obs])
})

test_that("m >= 2 returns a mimar_mi list of completed data.frames", {
  x <- mixed()
  out <- missknn(x, k = 3, m = 3, seed = 1, parallel_cores = 1L)
  expect_s3_class(out, "mimar_mi")
  expect_type(out, "list")
  expect_length(out, 3L)
  for (d in out) {
    expect_identical(class(d), "data.frame")
    expect_false(anyNA(d))
  }
  expect_false(isTRUE(all.equal(out[[1]], out[[2]])))
})

test_that("a model can be fitted on each imputed data set", {
  x <- mixed()
  fits <- lapply(missknn(x, k = 3, m = 2, seed = 1, parallel_cores = 1L),
                 function(d) stats::lm(a ~ b + g, data = d))
  expect_length(fits, 2L)
  expect_true(all(vapply(fits, function(f) all(is.finite(stats::coef(f))), logical(1))))
})

test_that("add_indicator appends one indicator per incomplete column", {
  x <- data.frame(a = c(1, NA, 3, 4), b = c(2, 4, 5, 6))
  out <- missknn(x, k = 1, add_indicator = TRUE)
  expect_identical(names(out), c("a", "b", "a_missing"))
  expect_identical(out$a_missing, c(0L, 1L, 0L, 0L))
})

test_that("seed makes multiple imputation reproducible", {
  x <- mixed()
  a <- missknn(x, k = 3, m = 2, seed = 7, parallel_cores = 1L)
  b <- missknn(x, k = 3, m = 2, seed = 7, parallel_cores = 1L)
  expect_identical(a, b)
})

test_that("invalid arguments are rejected", {
  x <- data.frame(a = c(1, NA, 3), b = c(2, 4, NA))
  expect_error(missknn(x, k = 0), "`k`")
  expect_error(missknn(x, m = 0), "`m`")
})

test_that("tuning attribute reports one row per variable", {
  x <- mixed()
  tu <- attr(missknn(x, k = 3), "tuning")
  expect_identical(tu$variable, names(x))
  expect_true(all(tu$k >= 1L))
  expect_true(all(tu$estimator[c(3, 4, 5)] == "vote"))
  expect_type(tu$global, "logical")
})

test_that("ordered factors and matrices are handled; tuning reports missing counts", {
  set.seed(3)
  d <- data.frame(
    x = c(1, NA, 3, 4, 5, 6, 7, NA),
    o = factor(c("lo", "hi", NA, "lo", "mid", "hi", "lo", "mid"),
               levels = c("lo", "mid", "hi"), ordered = TRUE)
  )
  out <- missknn(d, k = 2, seed = 1)
  expect_true(is.ordered(out$o))
  expect_identical(levels(out$o), c("lo", "mid", "hi"))
  expect_false(anyNA(out))
  expect_identical(class(out), "data.frame")
  expect_identical(attr(out, "tuning")$missing, c(2L, 1L))

  mx <- matrix(c(1, NA, 3, 4, 5, 6, 2, 3, NA, 5, 6, 7), ncol = 2)
  outm <- missknn(mx, k = 2, seed = 1)
  expect_identical(names(outm), c("V1", "V2"))
  expect_false(anyNA(outm))
})

test_that("blank strings are a category, not missing", {
  d <- data.frame(
    x = c(1, 2, NA, 4, 5, 6),
    s = c("u", "", NA, "u", "", "v"),
    stringsAsFactors = FALSE
  )
  out <- missknn(d, k = 2, seed = 1, add_indicator = TRUE)
  expect_identical(attr(out, "tuning")$missing, c(1L, 1L))
  expect_identical(out$s[c(1, 2, 4, 5, 6)], d$s[c(1, 2, 4, 5, 6)])
  expect_identical(out$s_missing, c(0L, 0L, 1L, 0L, 0L, 0L))
})
