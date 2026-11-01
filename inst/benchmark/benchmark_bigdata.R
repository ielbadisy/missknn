# Every method is timed on one thread: missRanger via num.threads = 1, missknn
# via RcppParallel's thread count (inherited by the per-run subprocesses).
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(ggplot2)
  library(knitr)
})

out_dir <- file.path("inst", "benchmark", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Each (n, method) combination is run in its own fresh `Rscript` subprocess.
# Running many large imputations back to back inside a single long-lived R
# session lets heap/GC growth from earlier iterations distort later timings
# (confirmed by profiling: near-zero self-time in missknn's own code even when
# wall time spikes) -- a fresh process per run avoids that confound and better
# reflects how a single big-data imputation call is actually used in practice.
runner_template <- "
suppressPackageStartupMessages({
  library(bench)
  library(missknn)
  library(missForest)
  library(missRanger)
  library(VIM)
})
set.seed(20260705)
# Equicorrelated Gaussian columns (pairwise correlation 0.5), so that
# neighbours carry information and accuracy is not capped at mean imputation.
simulate_data <- function(n, p = 5L, rho = 0.5) {
  z <- rnorm(n)
  x <- replicate(p, sqrt(rho) * z + sqrt(1 - rho) * rnorm(n), simplify = FALSE)
  names(x) <- paste0('x', seq_len(p))
  as.data.frame(x)
}
amputate <- function(df, prop = 0.2, seed) {
  set.seed(seed)
  out <- df
  n <- nrow(out)
  for (j in seq_len(ncol(out))) out[sample(n, floor(prop * n)), j] <- NA
  out
}
nrmse <- function(truth, imputed) {
  true_pooled <- unlist(truth)
  imp_pooled <- unlist(imputed)
  v <- stats::var(true_pooled)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  sqrt(mean((true_pooled - imp_pooled)^2) / v)
}

n <- %d
truth <- simulate_data(n)
miss <- amputate(truth, prop = 0.2, seed = 100L + %d)

method <- '%s'
call_fn <- switch(method,
  missknn = function() missknn(miss, k = 5L, m = 1L, seed = 1L),
  missForest = function() suppressWarnings(missForest::missForest(miss, verbose = FALSE)$ximp),
  missRanger = function() suppressWarnings(missRanger::missRanger(miss, verbose = 0, num.trees = 100, num.threads = 1)),
  `VIM::kNN` = function() suppressWarnings(as.data.frame(VIM::kNN(miss, k = 5L, imp_var = FALSE)))
)
# A single bench::mark() iteration is noisy for calls this fast: whether one
# GC pause happens to fall inside the timed call dominates the result. Average
# over several iterations for the cheap methods; the slow tree/distance-based
# methods are left at one iteration since repeating them would blow up total
# benchmark runtime for no diagnostic benefit.
reps <- if (method %%in%% c('missForest', 'missRanger', 'VIM::kNN')) 1L else 7L
timing <- bench::mark(call_fn(), iterations = reps, check = FALSE, filter_gc = FALSE)
out <- call_fn()
out <- as.data.frame(out)
err <- nrmse(truth[is.na(miss)], out[is.na(miss)])
cat(sprintf('%%.10f,%%.10f', as.numeric(timing$median), err))
"

# Every (n, method) cell is cached under output/bigdata_cells/, so the long
# runs can be resumed. Fast cells run one at a time; the slow cells
# (missForest and VIM::kNN beyond n = 5,000, hours each) then run side by side
# in separate processes (N_CORES, default 1).
cell_dir <- file.path(out_dir, "bigdata_cells")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
cores <- as.integer(Sys.getenv("N_CORES", "1"))

run_one <- function(n, i, method) {
  cache <- file.path(cell_dir, sprintf("%s_%d.csv", gsub("[^A-Za-z]", "", method), n))
  if (file.exists(cache)) {
    cached <- read.csv(cache, stringsAsFactors = FALSE)
    if (!"failed" %in% names(cached)) cached$failed <- FALSE
    if (!"failure" %in% names(cached)) cached$failure <- ""
    return(cached)
  }
  code <- sprintf(runner_template, n, i, method)
  script <- tempfile(fileext = ".R")
  writeLines(code, script)
  out <- tryCatch(system2("Rscript", script, stdout = TRUE, stderr = TRUE),
                  error = function(e) conditionMessage(e))
  parts <- suppressWarnings(tryCatch(as.numeric(strsplit(out[length(out)], ",")[[1]]),
                                     error = function(e) numeric()))
  failed <- length(parts) < 2L || any(!is.finite(parts[1:2])) ||
    !is.null(attr(out, "status"))
  res <- data.frame(n = n, method = method,
                    runtime_sec = if (failed) NA_real_ else parts[1],
                    nrmse = if (failed) NA_real_ else parts[2],
                    failed = failed,
                    failure = if (failed) paste(out, collapse = " | ") else "")
  write.csv(res, cache, row.names = FALSE)
  res
}

n_grid <- c(1000L, 5000L, 20000L, 100000L)
METHODS <- c("missknn", "missRanger", "missForest", "VIM::kNN")
cells <- expand.grid(i = seq_along(n_grid), method = METHODS, stringsAsFactors = FALSE)
cells$n <- n_grid[cells$i]
slow <- cells$method %in% c("missForest", "VIM::kNN") & cells$n > 5000L

run_cells <- function(k) run_one(cells$n[k], cells$i[k], cells$method[k])
results <- lapply(which(!slow), run_cells)
if (Sys.getenv("FAST_ONLY") == "1") quit(save = "no")  # cache the fast cells only
results <- c(results, parallel::mclapply(which(slow), run_cells, mc.cores = cores,
                                         mc.preschedule = FALSE))

bigdata_df <- do.call(rbind, results)
bigdata_df$method <- factor(bigdata_df$method, levels = c("missknn", "missForest", "missRanger", "VIM::kNN"))

runtime_plot <- ggplot(bigdata_df, aes(n, runtime_sec, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Runtime vs. Sample Size (log-log, one fresh process per run)",
    x = "n (rows)",
    y = "Seconds"
  ) +
  theme_minimal(base_size = 12)

nrmse_plot <- ggplot(bigdata_df, aes(n, nrmse, color = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Imputation NRMSE vs. Sample Size",
    x = "n (rows)",
    y = "NRMSE"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "bigdata_runtime_plot.png"), runtime_plot, width = 8, height = 5, dpi = 160)
ggsave(file.path(out_dir, "bigdata_mse_plot.png"), nrmse_plot, width = 8, height = 5, dpi = 160)
write.csv(bigdata_df, file.path(out_dir, "bigdata_results.csv"), row.names = FALSE)

table_md <- knitr::kable(
  transform(
    bigdata_df[, c("n", "method", "runtime_sec", "nrmse")],
    n = format(n, big.mark = ",", scientific = FALSE, trim = TRUE),
    runtime_sec = round(runtime_sec, 4), nrmse = round(nrmse, 5)
  ),
  format = "markdown",
  caption = "Runtime and NRMSE by method and sample size (one fresh Rscript process per run)"
)
failure_lines <- if (any(bigdata_df$failed)) {
  c("## Failed cells", "",
    sprintf("- n = %s, %s: %s", bigdata_df$n[bigdata_df$failed],
            bigdata_df$method[bigdata_df$failed], bigdata_df$failure[bigdata_df$failed]), "")
} else c("All method cells completed without recorded failures.", "")

report <- c(
  "# missknn Big-Data Benchmark Report",
  "",
  "Runtime and imputation accuracy (NRMSE) vs. sample size, from n = 1,000",
  "up to n = 100,000 (5 equicorrelated Gaussian columns, rho = 0.5, 20% MCAR per column). Each",
  "(n, method) run executes in its own fresh `Rscript` process, so timings",
  "are not distorted by heap growth accumulated across earlier, larger runs",
  "in a shared session. All four methods run across the full grid; the",
  "slow cells (`missForest` and `VIM::kNN` beyond n = 5,000) ran side by",
  "side in separate processes, one core each.",
  "",
  "## Figures",
  "",
  "![Runtime plot](bigdata_runtime_plot.png)",
  "",
  "![NRMSE plot](bigdata_mse_plot.png)",
  "",
  "## Results Table",
  "",
  table_md,
  "",
  failure_lines
)

writeLines(report, file.path(out_dir, "bigdata_benchmark_report.md"))
print(bigdata_df)
