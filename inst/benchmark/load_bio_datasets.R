## Fetch the three bioinformatics real datasets used in benchmark_bio.R,
## saved as .rds under inst/benchmark/data/. Mixed numeric + categorical
## tabular panels, sourced from genomics/proteomics/microarray repositories:
##   1. mice_protein   - UCI Mice Protein Expression (real MCAR-ish missingness
##                        from failed assay reads)
##   2. tcga_luad      - TCGA LUAD PanCancer Atlas clinical + genomic summary
##                        panel (mixed numeric/categorical), via cBioPortal
##   3. geo_lung       - GSE10072 lung cancer microarray (Affymetrix HG-U133A)
##                        expression + phenotype, via GEO series matrix
##
## The fetched .rds files are shipped alongside this script, so re-running it
## is only needed to verify provenance or refresh the data; benchmark_bio.R
## reads directly from inst/benchmark/data/ without needing network access.
## Run from the package root (Rscript inst/benchmark/load_bio_datasets.R).

suppressPackageStartupMessages({
  library(readxl)
  library(jsonlite)
  library(data.table)
})

out_dir <- file.path("inst", "benchmark", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cbio_api <- "https://www.cbioportal.org/api"

## ---- helper: paginated cBioPortal clinical-data fetch, wide-pivoted ----
## TCGA studies split clinical attributes across SAMPLE-level (e.g. mutation
## count, tumor purity) and PATIENT-level (e.g. age, sex, stage, survival)
## records; both are fetched and merged on patientId so the final panel
## isn't silently missing the demographic/outcome columns.
fetch_cbio_clinical <- function(study_id) {
  fetch_one <- function(level) {
    url <- sprintf("%s/studies/%s/clinical-data?clinicalDataType=%s&projection=SUMMARY&pageSize=20000",
                    cbio_api, study_id, level)
    dat <- as.data.frame(jsonlite::fromJSON(url, flatten = TRUE))
    if (!nrow(dat)) return(NULL)
    dat
  }
  samp <- fetch_one("SAMPLE")
  pat <- fetch_one("PATIENT")

  samp_map <- unique(samp[, c("sampleId", "patientId")])

  samp_wide <- as.data.frame(data.table::dcast(
    data.table::as.data.table(samp), patientId ~ clinicalAttributeId,
    value.var = "value", fun.aggregate = function(x) x[1]
  ))
  pat_wide <- as.data.frame(data.table::dcast(
    data.table::as.data.table(pat), patientId ~ clinicalAttributeId,
    value.var = "value", fun.aggregate = function(x) x[1]
  ))

  dup_cols <- intersect(setdiff(names(samp_wide), "patientId"), setdiff(names(pat_wide), "patientId"))
  pat_wide <- pat_wide[, setdiff(names(pat_wide), dup_cols)]

  wide <- merge(samp_map, samp_wide, by = "patientId")
  wide <- merge(wide, pat_wide, by = "patientId")
  wide$patientId <- NULL
  wide
}

## ---- 1. Mice Protein Expression (UCI ML repo #342) ----
message("[1/3] mice_protein: downloading from UCI ...")
mp_zip <- tempfile(fileext = ".zip")
download.file(
  "https://archive.ics.uci.edu/static/public/342/mice+protein+expression.zip",
  mp_zip, mode = "wb", quiet = TRUE
)
mp_dir <- tempfile()
unzip(mp_zip, exdir = mp_dir)
mp_xls <- list.files(mp_dir, pattern = "\\.xls$", full.names = TRUE)[1]
mice_protein <- as.data.frame(readxl::read_excel(mp_xls))
saveRDS(mice_protein, file.path(out_dir, "mice_protein.rds"))
message(sprintf("  -> mice_protein: n=%d p=%d, %.1f%% missing",
                 nrow(mice_protein), ncol(mice_protein),
                 100 * mean(is.na(mice_protein))))

## ---- 2. TCGA LUAD PanCancer Atlas clinical/genomic panel (cBioPortal) ----
message("[2/3] tcga_luad: fetching clinical/genomic panel from cBioPortal ...")
study3 <- "luad_tcga_pan_can_atlas_2018"
clin3 <- fetch_cbio_clinical(study3)
keep3 <- c("sampleId", "AGE", "SEX", "RACE", "AJCC_PATHOLOGIC_TUMOR_STAGE", "GRADE",
           "MUTATION_COUNT", "FRACTION_GENOME_ALTERED", "ANEUPLOIDY_SCORE",
           "TMB_NONSYNONYMOUS", "OS_MONTHS", "OS_STATUS")
tcga_luad <- clin3[, intersect(keep3, names(clin3))]
num_cols3 <- c("AGE", "MUTATION_COUNT", "FRACTION_GENOME_ALTERED", "ANEUPLOIDY_SCORE",
               "TMB_NONSYNONYMOUS", "OS_MONTHS")
for (cc in intersect(num_cols3, names(tcga_luad))) {
  tcga_luad[[cc]] <- suppressWarnings(as.numeric(tcga_luad[[cc]]))
}
fac_cols3 <- setdiff(names(tcga_luad), c("sampleId", num_cols3))
for (fc in fac_cols3) tcga_luad[[fc]] <- factor(tcga_luad[[fc]])

saveRDS(tcga_luad, file.path(out_dir, "tcga_luad.rds"))
message(sprintf("  -> tcga_luad: n=%d p=%d, %.1f%% missing",
                 nrow(tcga_luad), ncol(tcga_luad),
                 100 * mean(is.na(tcga_luad))))

## ---- 3. GSE10072 lung cancer microarray (GEO series matrix) ----
message("[3/3] geo_lung: downloading GSE10072 series matrix ...")
geo_gz <- tempfile(fileext = ".txt.gz")
## the NCBI FTP mirror rejects R's default download.file user-agent; curl with
## an explicit UA works
utils::download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10072/matrix/GSE10072_series_matrix.txt.gz",
  geo_gz, mode = "wb", quiet = TRUE, method = "curl",
  extra = '-A "Mozilla/5.0"'
)
geo_lines <- readLines(gzfile(geo_gz))

## phenotype rows are the "!Sample_*" header lines; expression matrix starts
## after "!series_matrix_table_begin"
pheno_idx <- grep("^!Sample_", geo_lines)
pheno_rows <- lapply(pheno_idx, function(i) {
  parts <- strsplit(geo_lines[i], "\t")[[1]]
  parts <- gsub('^"|"$', "", parts)
  parts
})
pheno_names <- vapply(pheno_rows, `[`, character(1), 1)
pheno_mat <- do.call(rbind, lapply(pheno_rows, `[`, -1))
rownames(pheno_mat) <- make.unique(sub("^!Sample_", "", pheno_names))
sample_ids <- pheno_mat["geo_accession", ]

## keep a few informative phenotype characteristics (cancer/normal status,
## gender, smoking) rather than the full free-text header
char_rows <- grep("characteristics_ch1", rownames(pheno_mat))
pheno_df <- data.frame(sampleId = as.character(sample_ids), stringsAsFactors = FALSE)
for (r in char_rows) {
  vals <- pheno_mat[r, ]
  tag <- sub(":.*$", "", vals[1])
  tag <- make.names(trimws(tag))
  pheno_df[[tag]] <- trimws(sub("^[^:]*:", "", vals))
}

begin_idx <- grep("^!series_matrix_table_begin", geo_lines)
end_idx <- grep("^!series_matrix_table_end", geo_lines)
expr_tab <- read.delim(
  text = geo_lines[(begin_idx + 1):(end_idx - 1)],
  header = TRUE, check.names = FALSE, quote = "\""
)
## restrict to the top 100 most-variable probes to keep this a tabular-scale
## panel (consistent with the mixed-tabular framing), transpose to samples x genes
probe_var <- apply(expr_tab[, -1], 1, stats::var, na.rm = TRUE)
top_idx <- order(probe_var, decreasing = TRUE)[seq_len(min(100, length(probe_var)))]
expr_top <- expr_tab[top_idx, ]
expr_num <- t(expr_top[, -1])
colnames(expr_num) <- expr_top[[1]]
expr_df <- data.frame(sampleId = rownames(expr_num), expr_num,
                       row.names = NULL, check.names = FALSE)

geo_lung <- merge(pheno_df, expr_df, by = "sampleId")
num_cols4 <- setdiff(names(expr_df), "sampleId")
for (cc in num_cols4) geo_lung[[cc]] <- suppressWarnings(as.numeric(geo_lung[[cc]]))
fac_cols4 <- setdiff(names(pheno_df), "sampleId")
for (fc in fac_cols4) geo_lung[[fc]] <- factor(geo_lung[[fc]])

saveRDS(geo_lung, file.path(out_dir, "geo_lung.rds"))
message(sprintf("  -> geo_lung: n=%d p=%d, %.1f%% missing",
                 nrow(geo_lung), ncol(geo_lung),
                 100 * mean(is.na(geo_lung))))

message("Done. Saved: mice_protein.rds, tcga_luad.rds, geo_lung.rds")
