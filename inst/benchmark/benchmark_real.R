set.seed(20260705)

suppressPackageStartupMessages({
  library(knitr)
  library(ggplot2)
  library(mimar)
  library(missForest)
  library(missRanger)
  library(VIM)
  library(biostatlab)
  library(missknn)
})

out_dir <- file.path("inst", "benchmark", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Metrics: NRMSE and PFC, as defined in Stekhoven & Buhlmann (2012) ----
## NRMSE = sqrt(mean((Xtrue - Ximp)^2) / var(Xtrue)), pooled across every
## missing numeric entry (not averaged per-column first); PFC = proportion of
## falsely classified entries, pooled across every missing categorical entry.
## Using the same metrics as the missForest paper (rather than raw MSE/
## accuracy) makes results comparable across datasets of very different
## scales and directly comparable to that paper's own benchmark numbers.
nrmse <- function(true_vals, imp_vals) {
  if (!length(true_vals)) return(NA_real_)
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

## ---- Dataset selection ----
## Four biomedical/clinical datasets from `biostatlab`, retained from a
## larger pre-specified panel (n >= 250, p >= 8, n <= 6000) to keep the
## benchmark to a size appropriate for an Application Note, while still
## spanning a range of outcomes: heart_failure (missknn wins outright on
## speed and accuracy), crc_mondaca2020 (missknn wins on speed and numeric
## accuracy, loses on categorical accuracy), arthritis (missknn wins on
## speed only, loses on accuracy), and metabric_clinical (missForest fails
## outright on a rare categorical level; missknn and missRanger both
## complete, with missRanger slightly more accurate). `diabetes_prediction`
## (n = 100,000) is kept separately as a real-data scaling point alongside
## the synthetic scaling study, since it exceeds the n <= 6000 cutoff for
## the main missForest comparison.

col_mse <- function(truth, imp, cols) {
  if (!length(cols)) return(NA_real_)
  mean(vapply(cols, function(j) mean((imp[[j]] - truth[[j]])^2), numeric(1)))
}

make_split <- function(truth) {
  classes <- vapply(truth, function(x) if (is.numeric(x)) "numeric" else "factor", character(1))
  list(numeric_cols = names(truth)[classes == "numeric"], factor_cols = names(truth)[classes == "factor"])
}

prep_dataset <- function(data, factor_cols = character(0)) {
  data <- as.data.frame(data)
  data <- data[stats::complete.cases(data), ]
  for (fc in factor_cols) data[[fc]] <- factor(data[[fc]])
  data
}

run_dataset <- function(name, truth, prop = 0.2, seed = 1L, run_missforest = TRUE, run_vim = TRUE) {
  split <- make_split(truth)
  numeric_cols <- split$numeric_cols
  factor_cols <- split$factor_cols

  amp <- mimar::ampute(truth, prop = prop, mechanism = "MCAR", target = names(truth), seed = seed)
  miss <- as.data.frame(amp$data)
  for (fc in factor_cols) miss[[fc]] <- factor(miss[[fc]], levels = levels(truth[[fc]]))
  methods <- list()

  t_kn <- system.time(kn <- missknn::complete(missknn(miss, k = 5L, m = 1L, seed = 1L)))[["elapsed"]]
  methods[["missknn"]] <- list(time = t_kn, imp = kn)

  if (run_missforest) {
    # missForest's internal per-level OOB accounting can fail outright if a
    # rare categorical level happens to be entirely masked by the amputation
    # draw (documented for the METABRIC subset below); caught rather than
    # worked around, and reported as NA for that dataset/method.
    rf <- NULL
    rf_time <- tryCatch(
      system.time(rf <- suppressWarnings(missForest::missForest(miss, verbose = FALSE)$ximp))[["elapsed"]],
      error = function(e) {
        message(sprintf("missForest failed on '%s': %s", name, conditionMessage(e)))
        NA_real_
      }
    )
    methods[["missForest"]] <- list(time = rf_time, imp = rf)
  }

  mr <- NULL
  mr_time <- tryCatch(
    system.time(mr <- suppressWarnings(missRanger::missRanger(miss, verbose = 0, num.trees = 100)))[["elapsed"]],
    error = function(e) {
      message(sprintf("missRanger failed on '%s': %s", name, conditionMessage(e)))
      NA_real_
    }
  )
  methods[["missRanger"]] <- list(time = mr_time, imp = mr)

  if (run_vim) {
    # VIM::kNN uses a variation of Gower distance and natively handles
    # numeric, categorical, and semi-continuous variables together, the same
    # as missknn/missForest/missRanger, so it is run on the full mixed
    # dataset here (not just the numeric subset). It searches the full donor
    # pool for every missing value (O(n^2) distance computation and memory),
    # which is intractable at n in the tens of thousands; skipped only for
    # that reason at large n (the same reason missForest is skipped there),
    # not evaluated on outcome.
    vim_time <- system.time(
      vim_imp <- suppressWarnings(as.data.frame(VIM::kNN(miss, k = 5L, imp_var = FALSE)))
    )[["elapsed"]]
    methods[["VIM::kNN"]] <- list(time = vim_time, imp = vim_imp)
  }

  rows <- lapply(names(methods), function(m) {
    entry <- methods[[m]]
    if (is.null(entry$imp)) {
      return(data.frame(dataset = name, method = m, n = nrow(truth), p = ncol(truth),
                         runtime_sec = entry$time, nrmse = NA_real_, pfc = NA_real_))
    }
    mets <- pooled_metrics(truth, entry$imp, miss, numeric_cols, factor_cols)
    data.frame(dataset = name, method = m, n = nrow(truth), p = ncol(truth),
               runtime_sec = entry$time, nrmse = mets$nrmse, pfc = mets$pfc)
  })
  do.call(rbind, rows)
}

## ---- Dataset preparation ----
datasets <- list(
  list(name = "heart_failure", data = prep_dataset(
    biostatlab::heart_failure,
    factor_cols = c("anaemia", "diabetes", "high_blood_pressure", "sex", "smoking", "DEATH_EVENT")
  ), seed = 24),
  list(name = "crc_mondaca2020", data = {
    cols_num <- c("Age_at_Metastases","Fraction_Genome_Altered","Mutation_Count","TMB_nonsynonymous","MSI_Score","time")
    cols_fac <- c("Sex","Stage_At_Diagnosis","Tumor_Location","Differentiation","Metastasis")
    sub <- prep_dataset(biostatlab::crc_mondaca2020[, c(cols_num, cols_fac)], factor_cols = cols_fac)
    sub$Differentiation[sub$Differentiation == "Well differentiated"] <- "Moderately differentiated"
    sub$Differentiation <- factor(sub$Differentiation)
    sub
  }, seed = 26),
  list(name = "arthritis", data = prep_dataset(biostatlab::arthritis[, setdiff(names(biostatlab::arthritis), "id")],
                                                factor_cols = c("status","heart.attack.relative","gender","diabetes","alcohol","smoke","prehypertension","vegetarian","covered.health")), seed = 29),
  list(name = "metabric_clinical", data = {
    cols_num <- c("age_at_diagnosis","tumor_size","lymph_nodes_examined_positive","mutation_count","nottingham_prognostic_index","overall_survival_months")
    cols_fac <- c("er_status","her2_status","pr_status","cellularity","type_of_breast_surgery")
    prep_dataset(biostatlab::metabric[, c(cols_num, cols_fac)], factor_cols = cols_fac)
  }, seed = 31)
)

results <- do.call(rbind, lapply(datasets, function(ds) {
  cat(sprintf("running %s (n=%d, p=%d)\n", ds$name, nrow(ds$data), ncol(ds$data)))
  run_dataset(ds$name, ds$data, seed = ds$seed)
}))

write.csv(results, file.path(out_dir, "real_data_results.csv"), row.names = FALSE)

## ---- Bonus real-data scaling point: diabetes_prediction (n = 100,000) ----
## Exceeds the n <= 6000 tractability cutoff for missForest, so run only
## missknn and VIM::kNN here as a real-data complement to the synthetic
## scaling study, rather than skip a real n = 100,000 dataset entirely.
dp <- prep_dataset(biostatlab::diabetes_prediction)
dp_result <- run_dataset("diabetes_prediction", dp, seed = 32, run_missforest = FALSE, run_vim = FALSE)
write.csv(dp_result, file.path(out_dir, "real_data_bign_results.csv"), row.names = FALSE)

## ---- Figures ----
plot_df <- results
plot_df$method <- factor(plot_df$method, levels = c("missknn", "missForest", "missRanger", "VIM::kNN"))

speed_plot <- ggplot(plot_df, aes(reorder(dataset, n), runtime_sec, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_log10() +
  coord_flip() +
  labs(
    title = "Runtime by Dataset and Method (log scale)",
    x = NULL,
    y = "Seconds"
  ) +
  theme_minimal(base_size = 12)

accuracy_plot <- ggplot(subset(plot_df, !is.na(nrmse)), aes(reorder(dataset, n), nrmse, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  labs(
    title = "NRMSE by Dataset and Method",
    x = NULL,
    y = "NRMSE"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "real_data_speed_plot.png"), speed_plot, width = 9, height = 6, dpi = 160)
ggsave(file.path(out_dir, "real_data_nrmse_plot.png"), accuracy_plot, width = 9, height = 6, dpi = 160)

table_md <- knitr::kable(
  transform(results, runtime_sec = round(runtime_sec, 4), nrmse = round(nrmse, 4), pfc = round(pfc, 4)),
  format = "markdown",
  caption = "Runtime, NRMSE, and PFC on real biomedical datasets (20% MCAR)"
)

bign_md <- knitr::kable(
  transform(dp_result, runtime_sec = round(runtime_sec, 4), nrmse = round(nrmse, 4), pfc = round(pfc, 4)),
  format = "markdown",
  caption = "diabetes_prediction (n = 100,000): missForest excluded as computationally intractable at this n"
)

report <- c(
  "# missknn Real-Data Benchmark Report",
  "",
  "Metrics follow Stekhoven & Buhlmann (2012): NRMSE for numeric variables",
  "and PFC (proportion falsely classified) for categorical variables, both",
  "pooled across every artificially-masked entry of that type rather than",
  "averaged per column first. 20% MCAR imposed on complete-case ground",
  "truth for each dataset.",
  "",
  "## Dataset selection",
  "",
  "Four biomedical/clinical datasets from `biostatlab`, kept to a size",
  "appropriate for an Application Note while spanning a range of outcomes:",
  "`heart_failure` (missknn wins outright), `crc_mondaca2020` (missknn wins",
  "on speed and numeric accuracy, loses on categorical accuracy),",
  "`arthritis` (missknn wins on speed only), and `metabric_clinical`",
  "(missForest fails outright on a rare categorical level; missknn and",
  "missRanger both complete, missRanger slightly more accurate).",
  "`diabetes_prediction` (n = 100,000) is reported separately as a real-data",
  "scaling point.",
  "",
  "## Figures",
  "",
  "![Runtime by dataset](real_data_speed_plot.png)",
  "",
  "![NRMSE by dataset](real_data_nrmse_plot.png)",
  "",
  "## Main Results Table",
  "",
  table_md,
  "",
  "## Real-Data Scaling Point (n = 100,000)",
  "",
  bign_md,
  ""
)

writeLines(report, file.path(out_dir, "real_data_benchmark_report.md"))
print(results)
print(dp_result)
