set.seed(20260802)

suppressPackageStartupMessages({
  library(knitr)
  library(ggplot2)
  library(mimar)
  library(missForest)
  library(missRanger)
  library(VIM)
  library(missknn)
})

data_dir <- file.path("inst", "benchmark", "data")
out_dir <- file.path("inst", "benchmark", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Metrics: NRMSE and PFC, as in Stekhoven & Buhlmann (2012) ----
## Same pooled-metric convention as the CIBM submission's benchmark_real.R:
## NRMSE/PFC pooled across every artificially-masked entry of that type
## across the whole dataset, not averaged per-column first.
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

run_dataset <- function(name, truth, prop = 0.3, seed = 1L) {
  split <- make_split(truth)
  numeric_cols <- split$numeric_cols
  factor_cols <- split$factor_cols

  amp <- mimar::ampute(truth, prop = prop, mechanism = "MCAR", target = names(truth), seed = seed)
  miss <- as.data.frame(amp$data)
  for (fc in factor_cols) miss[[fc]] <- factor(miss[[fc]], levels = levels(truth[[fc]]))
  methods <- list()

  t_kn <- system.time(kn <- missknn::complete(missknn(miss, k = 5L, m = 1L, seed = 1L)))[["elapsed"]]
  methods[["missknn"]] <- list(time = t_kn, imp = kn)

  rf <- NULL
  rf_time <- tryCatch(
    system.time(rf <- suppressWarnings(missForest::missForest(miss, verbose = FALSE)$ximp))[["elapsed"]],
    error = function(e) {
      message(sprintf("missForest failed on '%s': %s", name, conditionMessage(e)))
      NA_real_
    }
  )
  methods[["missForest"]] <- list(time = rf_time, imp = rf)

  mr <- NULL
  mr_time <- tryCatch(
    system.time(mr <- suppressWarnings(missRanger::missRanger(miss, verbose = 0, num.trees = 100, num.threads = 1)))[["elapsed"]],
    error = function(e) {
      message(sprintf("missRanger failed on '%s': %s", name, conditionMessage(e)))
      NA_real_
    }
  )
  methods[["missRanger"]] <- list(time = mr_time, imp = mr)

  vim_time <- tryCatch(
    system.time(vim_imp <- suppressWarnings(as.data.frame(VIM::kNN(miss, k = 5L, imp_var = FALSE))))[["elapsed"]],
    error = function(e) {
      message(sprintf("VIM::kNN failed on '%s': %s", name, conditionMessage(e)))
      NA_real_
    }
  )
  methods[["VIM::kNN"]] <- list(time = vim_time, imp = if (exists("vim_imp")) vim_imp else NULL)

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
## mice_protein / tcga_luad: standard complete-case ground truth, then 30%
## MCAR imposed (same convention as the CIBM real-data panel).
## geo_lung: top-100-variance probes already 0% missing as downloaded.

mp <- readRDS(file.path(data_dir, "mice_protein.rds"))
mp_factor_cols <- c("Genotype", "Treatment", "Behavior", "class")
mp <- prep_dataset(mp[, setdiff(names(mp), "MouseID")], factor_cols = mp_factor_cols)

tl <- readRDS(file.path(data_dir, "tcga_luad.rds"))
tl_factor_cols <- c("SEX", "RACE", "AJCC_PATHOLOGIC_TUMOR_STAGE", "OS_STATUS")
tl <- prep_dataset(tl[, setdiff(names(tl), "sampleId")], factor_cols = tl_factor_cols)

gl <- readRDS(file.path(data_dir, "geo_lung.rds"))
gl_factor_cols <- names(gl)[vapply(gl, is.factor, logical(1))]
gl <- prep_dataset(gl[, setdiff(names(gl), "sampleId")], factor_cols = gl_factor_cols)

datasets <- list(
  list(name = "mice_protein", data = mp, seed = 41),
  list(name = "tcga_luad", data = tl, seed = 47),
  list(name = "geo_lung", data = gl, seed = 53)
)

results <- do.call(rbind, lapply(datasets, function(ds) {
  cat(sprintf("running %s (n=%d, p=%d)\n", ds$name, nrow(ds$data), ncol(ds$data)))
  run_dataset(ds$name, ds$data, seed = ds$seed)
}))

write.csv(results, file.path(out_dir, "bio_data_results.csv"), row.names = FALSE)

## ---- Figures ----
plot_df <- results
plot_df$method <- factor(plot_df$method, levels = c("missknn", "missForest", "missRanger", "VIM::kNN"))

speed_plot <- ggplot(plot_df, aes(reorder(dataset, n), runtime_sec, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Runtime by Bioinformatics Dataset and Method (log scale)", x = NULL, y = "Seconds") +
  theme_minimal(base_size = 12)

accuracy_plot <- ggplot(subset(plot_df, !is.na(nrmse)), aes(reorder(dataset, n), nrmse, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  labs(title = "NRMSE by Bioinformatics Dataset and Method", x = NULL, y = "NRMSE") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "bio_data_speed_plot.png"), speed_plot, width = 9, height = 6, dpi = 160)
ggsave(file.path(out_dir, "bio_data_nrmse_plot.png"), accuracy_plot, width = 9, height = 6, dpi = 160)

table_md <- knitr::kable(
  transform(results, runtime_sec = round(runtime_sec, 4), nrmse = round(nrmse, 4), pfc = round(pfc, 4)),
  format = "markdown",
  caption = "Runtime, NRMSE, and PFC on bioinformatics datasets (30% MCAR)"
)

report <- c(
  "# missknn Bioinformatics Benchmark Report",
  "",
  "Metrics follow Stekhoven & Buhlmann (2012): NRMSE for numeric variables",
  "and PFC (proportion falsely classified) for categorical variables, both",
  "pooled across every artificially-masked entry of that type rather than",
  "averaged per column first. 30% MCAR imposed on complete-case/fully-observed",
  "ground truth for each dataset.",
  "",
  "## Dataset selection",
  "",
  "Three datasets sourced from genomics/proteomics/microarray repositories",
  "(not the `biostatlab` clinical panel used in the CIBM submission), to",
  "match a Bioinformatics-journal audience:",
  "",
  "- `mice_protein`: UCI Mice Protein Expression (protein/phosphoprotein",
  "  assay levels in mouse cortex; real assay-dropout missingness).",
  "- `tcga_luad`: TCGA Lung Adenocarcinoma PanCancer Atlas clinical/genomic",
  "  summary panel (age, stage, grade, mutation burden, survival), via",
  "  cBioPortal.",
  "- `geo_lung`: GSE10072 lung cancer microarray (Affymetrix HG-U133A),",
  "  top-100-variance probes + phenotype (cancer/normal, sex, smoking).",
  "",
  "## Figures",
  "",
  "![Runtime by dataset](bio_data_speed_plot.png)",
  "",
  "![NRMSE by dataset](bio_data_nrmse_plot.png)",
  "",
  "## Main Results Table",
  "",
  table_md,
  ""
)

writeLines(report, file.path(out_dir, "bio_data_benchmark_report.md"))
print(results)
