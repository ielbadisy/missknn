# Real-data scaling point: diabetes_prediction (n = 100,000), 30% MCAR, all
# four methods. Each method runs in its own fresh Rscript process and its
# result is cached under output/bign_cells/, so a long run can be resumed and
# the slow methods can run side by side (N_CORES, default 1).
# Every method is timed on one thread: missRanger via num.threads = 1, missknn
# via RcppParallel's thread count (inherited by the per-run subprocesses).
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressPackageStartupMessages(library(knitr))

out_dir <- file.path("inst", "benchmark", "output")
cell_dir <- file.path(out_dir, "bign_cells")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
cores <- as.integer(Sys.getenv("N_CORES", "1"))

METHODS <- c("missknn", "missForest", "missRanger", "VIM::kNN")

runner_template <- "
suppressPackageStartupMessages({
  library(missknn)
  library(missForest)
  library(missRanger)
  library(VIM)
})
nrmse <- function(t, i) { v <- stats::var(t); if (!is.finite(v) || v <= 0) NA_real_ else sqrt(mean((t - i)^2) / v) }
pfc <- function(t, i) if (!length(t)) NA_real_ else mean(as.character(t) != as.character(i))

truth <- as.data.frame(biostatlab::diabetes_prediction)
factor_cols <- c('gender', 'smoking_history', 'hypertension', 'heart_disease', 'diabetes')
for (fc in factor_cols) truth[[fc]] <- factor(truth[[fc]])
truth <- truth[stats::complete.cases(truth), ]
numeric_cols <- setdiff(names(truth), factor_cols)

amp <- mimar::ampute(truth, prop = 0.3, mechanism = 'MCAR', target = names(truth), seed = 32)
miss <- as.data.frame(amp$data)
for (fc in factor_cols) miss[[fc]] <- factor(miss[[fc]], levels = levels(truth[[fc]]))

method <- '%s'
t <- system.time(imp <- switch(method,
  missknn = missknn(miss, k = 5L, m = 1L, seed = 1L),
  missForest = suppressWarnings(missForest::missForest(miss, verbose = FALSE)$ximp),
  missRanger = suppressWarnings(missRanger::missRanger(miss, verbose = 0, num.trees = 100, num.threads = 1)),
  `VIM::kNN` = suppressWarnings(as.data.frame(VIM::kNN(miss, k = 5L, imp_var = FALSE)))
))[['elapsed']]
imp <- as.data.frame(imp)
nt <- unlist(lapply(numeric_cols, function(j) truth[[j]][is.na(miss[[j]])]))
ni <- unlist(lapply(numeric_cols, function(j) imp[[j]][is.na(miss[[j]])]))
ft <- unlist(lapply(factor_cols, function(j) as.character(truth[[j]])[is.na(miss[[j]])]))
fi <- unlist(lapply(factor_cols, function(j) as.character(imp[[j]])[is.na(miss[[j]])]))
cat(sprintf('%%d,%%d,%%.4f,%%.10f,%%.10f', nrow(truth), ncol(truth), t, nrmse(nt, ni), pfc(ft, fi)))
"

run_cell <- function(method) {
  cache <- file.path(cell_dir, paste0(gsub("[^A-Za-z]", "", method), ".csv"))
  if (file.exists(cache)) return(read.csv(cache))
  script <- tempfile(fileext = ".R")
  writeLines(sprintf(runner_template, method), script)
  out <- system2("Rscript", script, stdout = TRUE, stderr = FALSE)
  if (!length(out)) stop("run failed for ", method)
  parts <- as.numeric(strsplit(out[length(out)], ",")[[1]])
  res <- data.frame(dataset = "diabetes_prediction", method = method, n = parts[1], p = parts[2],
                    runtime_sec = parts[3], nrmse = parts[4], pfc = parts[5])
  write.csv(res, cache, row.names = FALSE)
  res
}

# The fast methods run alone so their timings are not disturbed; the slow
# ones then run side by side.
slow <- c("missForest", "VIM::kNN")
if (Sys.getenv("FAST_ONLY") == "1") {  # cache the fast cells only
  invisible(lapply(setdiff(METHODS, slow), run_cell))
  quit(save = "no")
}
dp_result <- do.call(rbind, c(lapply(setdiff(METHODS, slow), run_cell),
                              parallel::mclapply(slow, run_cell, mc.cores = cores,
                                                 mc.preschedule = FALSE)))
dp_result <- dp_result[match(METHODS, dp_result$method), ]
write.csv(dp_result, file.path(out_dir, "real_data_bign_results.csv"), row.names = FALSE)

report <- c(
  "# missknn Real-Data Scaling Point",
  "",
  "`diabetes_prediction` (n = 100,000, p = 9; gender, smoking history and the",
  "three binary indicators treated as categorical), 30% MCAR, all four",
  "methods, each in its own fresh Rscript process.",
  "",
  knitr::kable(transform(dp_result, runtime_sec = round(runtime_sec, 1), nrmse = round(nrmse, 4),
                         pfc = round(pfc, 4)), format = "markdown"),
  ""
)
writeLines(report, file.path(out_dir, "real_data_bign_report.md"))
print(dp_result)
