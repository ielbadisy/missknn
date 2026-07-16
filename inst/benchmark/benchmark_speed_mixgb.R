## Execution-time comparison: missknn vs. missForest, missRanger, mixgb.
##
## Single MCAR imputation pass at increasing n, measured with bench::mark
## (one iteration each, since missForest/mixgb are too slow to repeat many
## times at the largest n). Output: runtime table, CSV, and plot under
## inst/benchmark/output/.

set.seed(20260716)

suppressPackageStartupMessages({
  library(bench)
  library(ggplot2)
  library(knitr)
  library(mimar)
  library(missForest)
  library(missRanger)
  library(mixgb)
  library(missknn)
})

# All methods pinned to a single core/thread for a fair, single-imputation
# runtime comparison: missForest is serial by default (parallelize = "no"),
# missRanger/mixgb take explicit thread args, and missknn's RcppParallel
# layer (per-column k/estimator tuning + neighbor search) is capped via
# RCPP_PARALLEL_NUM_THREADS.
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
data.table::setDTthreads(1L)

out_dir <- file.path("inst", "benchmark", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

METHODS <- c("missknn", "missForest", "missRanger", "mixgb")

simulate_data <- function(n, p = 8L) {
  x <- replicate(p, rnorm(n), simplify = FALSE)
  names(x) <- paste0("x", seq_len(p))
  as.data.frame(x)
}

ampute <- function(df, prop = 0.2, seed) {
  amp <- suppressWarnings(
    mimar::ampute(df, prop = prop, mechanism = "MCAR", target = names(df), seed = seed)
  )
  as.data.frame(amp$data)
}

time_one <- function(miss, seed) {
  bench::mark(
    missknn = complete(missknn(miss, k = 5L, m = 1L, seed = seed, parallel_cores = 1L)),
    missForest = suppressWarnings(missForest::missForest(miss, verbose = FALSE, parallelize = "no")$ximp),
    missRanger = suppressWarnings(missRanger::missRanger(miss, verbose = 0, num.trees = 100, num.threads = 1)),
    mixgb = suppressMessages(suppressWarnings(
      mixgb::mixgb(miss, m = 1L, verbose = FALSE, xgb.params = list(nthread = 1))[[1]]
    )),
    iterations = 1L, check = FALSE
  )
}

sizes <- c(500L, 2000L, 10000L, 50000L)

results <- do.call(rbind, lapply(sizes, function(n) {
  df <- simulate_data(n)
  miss <- ampute(df, seed = n)
  timing <- time_one(miss, seed = n)
  data.frame(
    n = n,
    method = as.character(timing$expression),
    runtime_sec = as.numeric(timing$median)
  )
}))

write.csv(results, file.path(out_dir, "speed_mixgb_results.csv"), row.names = FALSE)

wide <- reshape(results, idvar = "n", timevar = "method", direction = "wide")
names(wide) <- sub("runtime_sec\\.", "", names(wide))
wide[METHODS] <- lapply(wide[METHODS], round, digits = 4)

p <- ggplot(results, aes(x = n, y = runtime_sec, color = method)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "n (rows)", y = "median runtime (sec, log scale)",
    title = "missknn vs. missForest, missRanger, mixgb: single-imputation runtime"
  ) +
  theme_minimal()
ggsave(file.path(out_dir, "speed_mixgb_plot.png"), p, width = 7, height = 4.5, dpi = 150)

report <- c(
  "# missknn vs. missForest, missRanger, mixgb: Execution Time",
  "",
  "Single MCAR imputation (20% per column), p = 8 numeric predictors,",
  "n scaling from 500 to 50,000. One measured run per (method, n) via",
  "`bench::mark(iterations = 1)`, since missForest/mixgb are too slow to",
  "repeat many times at the largest n.",
  "",
  "![Runtime plot](speed_mixgb_plot.png)",
  "",
  "## Runtime Table (seconds, median)",
  "",
  knitr::kable(wide, format = "markdown", row.names = FALSE),
  ""
)
writeLines(report, file.path(out_dir, "speed_mixgb_report.md"))
print(wide)
