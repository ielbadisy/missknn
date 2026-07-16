set.seed(20260705)

suppressPackageStartupMessages({
  library(bench)
  library(ggplot2)
  library(knitr)
  library(mimar)
  library(missForest)
  library(missRanger)
  library(VIM)
  library(missknn)
})

out_dir <- file.path("inst", "benchmark", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

METHODS <- c("missknn", "missForest", "missRanger", "VIM::kNN")

simulate_data <- function(n = 200L, p = 5L) {
  x <- replicate(p, rnorm(n), simplify = FALSE)
  names(x) <- paste0("x", seq_len(p))
  x <- as.data.frame(x)
  beta <- c(0.6, -0.4, 0.35, 0.25, -0.2)
  eta <- 0.5 + as.matrix(x) %*% beta
  y <- as.numeric(eta + rnorm(n, sd = 0.75))
  list(data = data.frame(y = y, x, check.names = FALSE), beta = beta)
}

ampute_predictors <- function(df, prop = 0.2, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  x <- df[, setdiff(names(df), "y"), drop = FALSE]
  amp <- suppressWarnings(mimar::ampute(x, prop = prop, mechanism = "MCAR", target = names(x), seed = seed))
  list(
    complete = data.frame(y = df$y, amp$data, check.names = FALSE),
    truth = df
  )
}

fit_lm <- function(df) {
  stats::coef(stats::lm(y ~ ., data = df))
}

# NRMSE pooled across every imputed (originally-missing) predictor entry, as
# in Stekhoven & Buhlmann (2012); consistent with the metric used throughout
# inst/benchmark/benchmark_real.R.
nrmse <- function(true_vals, imp_vals) {
  v <- stats::var(true_vals)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  sqrt(mean((true_vals - imp_vals)^2) / v)
}

impute_methods <- function(df_miss, truth_x, seed) {
  x_miss <- df_miss[, setdiff(names(df_miss), "y"), drop = FALSE]

  timings <- bench::mark(
    missknn = complete(missknn(x_miss, k = 5L, m = 1L, seed = seed)),
    missForest = suppressWarnings(missForest::missForest(x_miss, verbose = FALSE)$ximp),
    missRanger = suppressWarnings(missRanger::missRanger(x_miss, verbose = 0, num.trees = 100, num.threads = 1)),
    `VIM::kNN` = suppressWarnings(as.data.frame(VIM::kNN(x_miss, k = 5L, imp_var = FALSE))),
    iterations = 1,
    check = FALSE
  )

  completed <- list(
    missknn = complete(missknn(x_miss, k = 5L, m = 1L, seed = seed)),
    missForest = suppressWarnings(missForest::missForest(x_miss, verbose = FALSE)$ximp),
    missRanger = suppressWarnings(missRanger::missRanger(x_miss, verbose = 0, num.trees = 100, num.threads = 1)),
    `VIM::kNN` = suppressWarnings(as.data.frame(VIM::kNN(x_miss, k = 5L, imp_var = FALSE)))
  )

  na_mask <- as.data.frame(is.na(x_miss))
  nrmse_by_method <- vapply(METHODS, function(m) {
    imp <- completed[[m]]
    true_pooled <- unlist(lapply(names(truth_x), function(j) truth_x[[j]][na_mask[[j]]]))
    imp_pooled <- unlist(lapply(names(truth_x), function(j) imp[[j]][na_mask[[j]]]))
    nrmse(true_pooled, imp_pooled)
  }, numeric(1))

  list(timings = timings, completed = completed, nrmse = nrmse_by_method)
}

score_methods <- function(beta_full, beta_hat) {
  common <- intersect(names(beta_full), names(beta_hat))
  diff <- beta_hat[common] - beta_full[common]
  data.frame(
    term = common,
    bias = diff,
    mse = diff^2,
    row.names = NULL
  )
}

results <- list()
coef_results <- list()
runtime_results <- list()
nrmse_results <- list()

for (s in seq_len(30L)) {
  sim <- simulate_data()
  amp <- ampute_predictors(sim$data, prop = 0.2, seed = s)
  full_beta <- fit_lm(amp$truth)
  truth_x <- sim$data[, setdiff(names(sim$data), "y"), drop = FALSE]

  imp <- impute_methods(amp$complete, truth_x, seed = s)
  timing_tbl <- as.data.frame(imp$timings)
  timing_tbl$method <- METHODS
  timing_tbl$simulation <- s
  runtime_results[[s]] <- timing_tbl[, c("simulation", "method", "median", "mem_alloc")]

  nrmse_results[[s]] <- data.frame(simulation = s, method = METHODS, nrmse = imp$nrmse, row.names = NULL)

  for (method in METHODS) {
    completed <- data.frame(y = amp$complete$y, imp$completed[[method]], check.names = FALSE)
    beta_hat <- fit_lm(completed)
    sc <- score_methods(full_beta, beta_hat)
    sc$simulation <- s
    sc$method <- method
    coef_results[[length(coef_results) + 1L]] <- sc
  }
}

runtime_df <- do.call(rbind, runtime_results)
runtime_df$runtime_sec <- as.numeric(runtime_df$median)
runtime_df$method <- factor(runtime_df$method, levels = METHODS)

coef_df <- do.call(rbind, coef_results)
coef_df$method <- factor(coef_df$method, levels = METHODS)

nrmse_df <- do.call(rbind, nrmse_results)
nrmse_df$method <- factor(nrmse_df$method, levels = METHODS)

runtime_summary <- aggregate(runtime_sec ~ method, runtime_df, mean)
runtime_summary$sd_runtime <- aggregate(runtime_sec ~ method, runtime_df, sd)$runtime_sec

coef_summary <- aggregate(cbind(abs_bias = abs(bias), mse = mse) ~ method, coef_df, mean)

coef_term_summary <- aggregate(cbind(abs_bias = abs(bias), mse = mse) ~ method + term, coef_df, mean)

nrmse_summary <- aggregate(nrmse ~ method, nrmse_df, mean)

runtime_plot <- ggplot(runtime_df, aes(method, runtime_sec, fill = method)) +
  geom_boxplot(width = 0.65, outlier.alpha = 0.35) +
  scale_y_log10() +
  labs(
    title = "Benchmark Runtime Across 30 Simulations (log scale)",
    x = NULL,
    y = "Seconds"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

bias_plot <- ggplot(coef_term_summary, aes(term, abs_bias, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Mean Absolute Coefficient Bias by Method",
    x = "Coefficient",
    y = "Absolute bias"
  ) +
  theme_minimal(base_size = 12)

mse_plot <- ggplot(coef_summary, aes(method, mse, fill = method)) +
  geom_col(width = 0.7) +
  labs(
    title = "Mean Coefficient MSE by Method",
    x = NULL,
    y = "MSE"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

nrmse_plot <- ggplot(nrmse_summary, aes(method, nrmse, fill = method)) +
  geom_col(width = 0.7) +
  labs(
    title = "Mean NRMSE by Method",
    x = NULL,
    y = "NRMSE"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "runtime_plot.png"), runtime_plot, width = 8, height = 5, dpi = 160)
ggsave(file.path(out_dir, "bias_plot.png"), bias_plot, width = 9, height = 5, dpi = 160)
ggsave(file.path(out_dir, "mse_plot.png"), mse_plot, width = 8, height = 5, dpi = 160)
ggsave(file.path(out_dir, "nrmse_plot.png"), nrmse_plot, width = 8, height = 5, dpi = 160)

write.csv(runtime_df, file.path(out_dir, "runtime_results.csv"), row.names = FALSE)
write.csv(coef_df, file.path(out_dir, "coefficient_results.csv"), row.names = FALSE)
write.csv(nrmse_df, file.path(out_dir, "nrmse_results.csv"), row.names = FALSE)

runtime_md <- knitr::kable(
  transform(runtime_summary, runtime_sec = round(runtime_sec, 4), sd_runtime = round(sd_runtime, 4)),
  format = "markdown",
  caption = "Mean runtime across 30 simulations"
)

nrmse_md <- knitr::kable(
  transform(nrmse_summary, nrmse = round(nrmse, 4)),
  format = "markdown",
  caption = "Mean NRMSE (imputed predictors, pooled) across 30 simulations"
)

coef_md <- knitr::kable(
  transform(coef_summary, abs_bias = round(abs_bias, 4), mse = round(mse, 4)),
  format = "markdown",
  caption = "Mean downstream OLS coefficient error across 30 simulations"
)

term_md <- knitr::kable(
  transform(coef_term_summary, abs_bias = round(abs_bias, 4), mse = round(mse, 4)),
  format = "markdown",
  caption = "Coefficient-level absolute bias and MSE"
)

report <- c(
  "# missknn Benchmark Report",
  "",
  "30 simulations, MCAR amputation from `mimar`, `bench::mark` for runtime.",
  "Two kinds of accuracy are reported: NRMSE (Stekhoven & Buhlmann 2012),",
  "measuring imputation accuracy on the predictors directly, and downstream",
  "OLS coefficient bias/MSE, measuring how imputation error propagates into",
  "a regression fit on the imputed data -- the original evaluation design",
  "from the missForest paper.",
  "",
  "## Figures",
  "",
  "![Runtime plot](runtime_plot.png)",
  "",
  "![NRMSE plot](nrmse_plot.png)",
  "",
  "![Bias plot](bias_plot.png)",
  "",
  "![MSE plot](mse_plot.png)",
  "",
  "## Runtime Table",
  "",
  runtime_md,
  "",
  "## NRMSE Table",
  "",
  nrmse_md,
  "",
  "## Downstream Coefficient Error Table",
  "",
  coef_md,
  "",
  "## Coefficient-Level Table",
  "",
  term_md,
  ""
)

writeLines(report, file.path(out_dir, "benchmark_report.md"))

print(runtime_summary)
print(nrmse_summary)
print(coef_summary)
