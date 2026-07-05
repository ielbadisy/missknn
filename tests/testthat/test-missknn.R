test_that("single imputation returns a completed object", {
  x <- data.frame(a = c(1, NA, 3), b = c(2, 4, NA))
  imp <- missknn::missknn(x, k = 1)
  out <- missknn::complete(imp)
  expect_true(is.data.frame(out) || data.table::is.data.table(out))
  expect_false(anyNA(out))
})

test_that("multiple imputation returns a list", {
  x <- data.frame(a = c(1, NA, 3), b = c(2, 4, NA))
  imp <- missknn::missknn(x, k = 1, m = 2, seed = 1, parallel_cores = 1L)
  out <- missknn::complete(imp)
  expect_type(out, "list")
  expect_equal(length(out), 2L)
})
