set.seed(20260705)

suppressPackageStartupMessages({
  library(knitr)
  library(missForest)
  library(missRanger)
  library(VIM)
  library(mimar)
  library(missknn)
})

out_dir <- file.path("inst", "benchmark", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

METHODS <- c("missknn", "missForest", "missRanger", "VIM::kNN")

nrmse <- function(true_vals, imp_vals) {
  v <- stats::var(true_vals)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  sqrt(mean((true_vals - imp_vals)^2) / v)
}
pfc <- function(true_vals, imp_vals) {
  if (!length(true_vals)) return(NA_real_)
  mean(as.character(true_vals) != as.character(imp_vals))
}
pooled_metrics <- function(truth, imp, miss, numeric_cols, factor_cols) {
  num_true <- unlist(lapply(numeric_cols, function(j) truth[[j]][is.na(miss[[j]])]))
  num_imp <- unlist(lapply(numeric_cols, function(j) imp[[j]][is.na(miss[[j]])]))
  fac_true <- unlist(lapply(factor_cols, function(j) as.character(truth[[j]])[is.na(miss[[j]])]))
  fac_imp <- unlist(lapply(factor_cols, function(j) as.character(imp[[j]])[is.na(miss[[j]])]))
  list(nrmse = nrmse(num_true, num_imp), pfc = pfc(fac_true, fac_imp))
}

run_all_methods <- function(miss, numeric_cols, factor_cols) {
  list(
    missknn = missknn::complete(missknn(miss, k = 5L, m = 1L, seed = 1L)),
    missForest = suppressWarnings(missForest::missForest(miss, verbose = FALSE)$ximp),
    missRanger = suppressWarnings(missRanger::missRanger(miss, verbose = 0, num.trees = 100)),
    `VIM::kNN` = suppressWarnings(as.data.frame(VIM::kNN(miss, k = 5L, imp_var = FALSE)))
  )
}

## ---- Study 1: correlated numeric data -------------------------------------
## Checks whether missknn's local-regression / global-fallback tuning stays
## competitive once predictors carry genuine linear correlation, rather than
## only on the independent-noise case.
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
imps_cor <- run_all_methods(miss_cor, names(truth_cor), character(0))

correlated_results <- do.call(rbind, lapply(METHODS, function(m) {
  mets <- pooled_metrics(truth_cor, imps_cor[[m]], miss_cor, names(truth_cor), character(0))
  data.frame(method = m, nrmse = mets$nrmse)
}))

## ---- Study 2: mixed numeric + categorical data ----------------------------
n <- 400
x1 <- rnorm(n)
x2 <- 0.7 * x1 + rnorm(n, sd = 0.7)
x3 <- rnorm(n)
g <- factor(ifelse(x1 + rnorm(n, sd = 0.5) > 0, "hi", "lo"))
truth_mix <- data.frame(x1 = x1, x2 = x2, x3 = x3, g = g)

miss_mix <- truth_mix
for (j in c("x1", "x2", "x3")) miss_mix[sample(n, floor(0.2 * n)), j] <- NA
miss_mix$g[sample(n, floor(0.2 * n))] <- NA
imps_mix <- run_all_methods(miss_mix, c("x1", "x2", "x3"), "g")

mixed_results <- do.call(rbind, lapply(METHODS, function(m) {
  mets <- pooled_metrics(truth_mix, imps_mix[[m]], miss_mix, c("x1", "x2", "x3"), "g")
  data.frame(method = m, nrmse = mets$nrmse, pfc = mets$pfc)
}))

write.csv(correlated_results, file.path(out_dir, "mixed_correlated_results.csv"), row.names = FALSE)
write.csv(mixed_results, file.path(out_dir, "mixed_categorical_results.csv"), row.names = FALSE)

report <- c(
  "# missknn Correlated- and Mixed-Data Validation",
  "",
  "Metrics: NRMSE (numeric) and PFC (categorical), as in",
  "Stekhoven & Buhlmann (2012); see inst/benchmark/benchmark_real.R.",
  "",
  "## Study 1: correlated numeric data",
  "",
  "n = 300, p = 6 numeric variables drawn from an equicorrelated Gaussian",
  "(rho = 0.6), 20% MCAR per column.",
  "",
  knitr::kable(transform(correlated_results, nrmse = round(nrmse, 5)), format = "markdown"),
  "",
  "## Study 2: mixed numeric + categorical data",
  "",
  "n = 400, 3 correlated numeric variables plus 1 binary factor driven by one",
  "of them, 20% MCAR per column.",
  "",
  knitr::kable(
    transform(mixed_results, nrmse = round(nrmse, 5), pfc = round(pfc, 4)),
    format = "markdown"
  ),
  ""
)

writeLines(report, file.path(out_dir, "mixed_validation_report.md"))
print(correlated_results)
print(mixed_results)
