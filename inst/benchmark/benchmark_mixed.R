set.seed(20260705)

suppressPackageStartupMessages({
  library(knitr)
  library(missForest)
  library(mimar)
  library(missknn)
})

out_dir <- file.path("inst", "benchmark", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

col_mse <- function(truth, imp, cols) {
  mean(vapply(cols, function(j) mean((imp[[j]] - truth[[j]])^2), numeric(1)))
}

## ---- Study 1: correlated numeric data -------------------------------------
## Checks whether missknn's local-regression / global-fallback tuning stays
## competitive with missForest once predictors carry genuine linear
## correlation, rather than only on the independent-noise case.
simulate_correlated <- function(n = 300, p = 6, rho = 0.6) {
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  L <- chol(Sigma)
  z <- matrix(rnorm(n * p), n, p) %*% L
  colnames(z) <- paste0("x", seq_len(p))
  as.data.frame(z)
}

truth_cor <- simulate_correlated()
amp_cor <- mimar::ampute(truth_cor, prop = 0.2, mechanism = "MCAR", target = names(truth_cor), seed = 7)
miss_cor <- as.data.frame(amp_cor$data)

kn_cor <- missknn::complete(missknn(miss_cor, k = 5L, m = 1L, seed = 1L))
rf_cor <- suppressWarnings(missForest::missForest(miss_cor, verbose = FALSE)$ximp)

correlated_results <- data.frame(
  method = c("missknn", "missForest"),
  mse = c(
    col_mse(truth_cor, kn_cor, names(truth_cor)),
    col_mse(truth_cor, rf_cor, names(truth_cor))
  )
)

## ---- Study 2: mixed numeric + categorical data ----------------------------
## missknn and missForest both handle numeric and categorical targets from
## the same call.
n <- 400
x1 <- rnorm(n)
x2 <- 0.7 * x1 + rnorm(n, sd = 0.7)
x3 <- rnorm(n)
g <- factor(ifelse(x1 + rnorm(n, sd = 0.5) > 0, "hi", "lo"))
truth_mix <- data.frame(x1 = x1, x2 = x2, x3 = x3, g = g)

miss_mix <- truth_mix
for (j in c("x1", "x2", "x3")) miss_mix[sample(n, floor(0.2 * n)), j] <- NA
miss_mix$g[sample(n, floor(0.2 * n))] <- NA

kn_mix <- missknn::complete(missknn(miss_mix, k = 5L, m = 1L, seed = 1L))
rf_mix <- suppressWarnings(missForest::missForest(miss_mix, verbose = FALSE)$ximp)

mixed_results <- data.frame(
  method = c("missknn", "missForest"),
  numeric_mse = c(
    col_mse(truth_mix, kn_mix, c("x1", "x2", "x3")),
    col_mse(truth_mix, rf_mix, c("x1", "x2", "x3"))
  ),
  factor_accuracy = c(
    mean(kn_mix$g == truth_mix$g),
    mean(rf_mix$g == truth_mix$g)
  )
)

write.csv(correlated_results, file.path(out_dir, "mixed_correlated_results.csv"), row.names = FALSE)
write.csv(mixed_results, file.path(out_dir, "mixed_categorical_results.csv"), row.names = FALSE)

report <- c(
  "# missknn Correlated- and Mixed-Data Validation",
  "",
  "## Study 1: correlated numeric data",
  "",
  "n = 300, p = 6 numeric variables drawn from an equicorrelated Gaussian",
  "(rho = 0.6), 20% MCAR per column.",
  "",
  knitr::kable(transform(correlated_results, mse = round(mse, 5)), format = "markdown"),
  "",
  "## Study 2: mixed numeric + categorical data",
  "",
  "n = 400, 3 correlated numeric variables plus 1 binary factor driven by one",
  "of them, 20% MCAR per column. `missknn` and `missForest` handle numeric and",
  "categorical targets from the same call.",
  "",
  knitr::kable(
    transform(mixed_results, numeric_mse = round(numeric_mse, 5), factor_accuracy = round(factor_accuracy, 4)),
    format = "markdown"
  ),
  ""
)

writeLines(report, file.path(out_dir, "mixed_validation_report.md"))
print(correlated_results)
print(mixed_results)
