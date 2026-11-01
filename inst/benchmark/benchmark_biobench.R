# Secondary BioBench comparison. Run from the package root with
# Rscript inst/benchmark/benchmark_biobench.R. Results checkpoint per cell.
Sys.setenv(RCPP_PARALLEL_NUM_THREADS="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
suppressPackageStartupMessages({library(biobench); library(missknn); library(missForest); library(missRanger); library(VIM); library(R.utils)})
out <- "inst/benchmark/output/biobench"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
methods <- c("missknn","missForest","missRanger","VIM::kNN")
idx <- biobench::biobench_index(outcome="regression", missing=FALSE)
path <- file.path(out,"biobench_all_results.csv")
res <- if(file.exists(path)) read.csv(path,stringsAsFactors=FALSE) else data.frame()
if(nrow(res) && !"attempted" %in% names(res)) res$attempted <- TRUE
metric <- function(truth,imp,mask) {
  ni <- which(vapply(truth,is.numeric,logical(1))); ci <- setdiff(seq_along(truth),ni)
  a <- unlist(lapply(ni,function(j) truth[[j]][mask[,j]])); b <- unlist(lapply(ni,function(j) imp[[j]][mask[,j]]))
  nrmse <- if(length(a)>1 && var(a)>0) sqrt(mean((a-b)^2)/var(a)) else NA_real_
  a <- unlist(lapply(ci,function(j) as.character(truth[[j]][mask[,j]]))); b <- unlist(lapply(ci,function(j) as.character(imp[[j]][mask[,j]])))
  c(nrmse=nrmse,pfc=if(length(a)) mean(a!=b) else NA_real_)
}
fit <- function(m,x) switch(m,
  missknn=missknn::missknn(x,k=5L,m=1L,seed=1L),
  missForest=missForest::missForest(x,verbose=FALSE)$ximp,
  missRanger=missRanger::missRanger(x,verbose=0,num.trees=100,num.threads=1),
  `VIM::kNN`=VIM::kNN(x,k=5L,imp_var=FALSE))
for(di in seq_len(nrow(idx))) {
  nm <- idx$dataset[di]; e <- new.env(parent=emptyenv()); utils::data(list=nm,package="biobench",envir=e)
  dat <- as.data.frame(get(nm,envir=e)); dat <- dat[complete.cases(dat),,drop=FALSE]
  id <- vapply(dat,function(z) length(unique(z))==nrow(dat),logical(1)); dat <- dat[,!id,drop=FALSE]
  for(j in seq_along(dat)) if(is.character(dat[[j]])||is.logical(dat[[j]])) dat[[j]] <- factor(dat[[j]])
  # Include every non-identifier column in the known-beta regression. Factors
  # enter as indicator columns; all design columns use fixed complete-data
  # centering/scaling so the coefficient target is identical after imputation.
  design <- function(x) model.matrix(~ .,data=x)[,-1,drop=FALSE]
  xraw <- design(dat); xcenter <- colMeans(xraw); xscale <- apply(xraw,2,sd)
  xscale[!is.finite(xscale)|xscale==0] <- 1
  xx <- scale(xraw,center=xcenter,scale=xscale); xx[!is.finite(xx)] <- 0
  beta <- seq(.4,-.2,length.out=ncol(xx))
  set.seed(8000+di); y <- drop(xx%*%beta)+rnorm(nrow(xx))
  for(rep in 1:30) {
    seed <- 20261101L+di*10000L+rep; set.seed(seed)
    truth <- dat; miss <- truth
    mask <- matrix(runif(nrow(truth) * ncol(truth)) < .30,
      nrow = nrow(truth), ncol = ncol(truth))
    stopifnot(identical(dim(mask), dim(as.matrix(truth))), all(colSums(mask) > 0L))
    for(j in seq_along(truth)) miss[[j]][mask[,j]] <- NA
    for(m in methods) {
      key <- paste(nm,rep,m,sep="|"); if(nrow(res)&&key%in%paste(res$dataset,res$replicate,res$method,sep="|")) next
      msg <- ""; tic <- proc.time()[["elapsed"]]
      attempted <- TRUE
      agefat_timeout <- nm == "agefat" && m == "missForest" && rep == 12L &&
        any(res$dataset == "agefat" & res$method == "missForest" &
              res$replicate == 12L & res$failed)
      if(nm == "agefat" && m == "missForest" && rep > 12L && agefat_timeout) {
        attempted <- FALSE
        msg <- "not attempted after the replicate-12 missForest failure on agefat"
        imp <- NULL
      }
      else imp <- tryCatch(R.utils::withTimeout(fit(m,miss),timeout=120,onTimeout="error"),error=function(e){msg<<-conditionMessage(e);NULL})
      elapsed <- proc.time()[["elapsed"]]-tic; met <- if(is.null(imp)) c(nrmse=NA,pfc=NA) else metric(truth,as.data.frame(imp),mask)
      bias <- NA_real_
      if(!is.null(imp)) { zraw <- design(as.data.frame(imp)); zraw <- zraw[,colnames(xraw),drop=FALSE]
        z <- scale(zraw,center=xcenter,scale=xscale); z[!is.finite(z)] <- 0
        b <- tryCatch(coef(lm(y~z))[-1],error=function(e) rep(NA_real_,length(beta))); bias <- mean(abs(b-beta),na.rm=TRUE) }
      res <- rbind(res,data.frame(dataset=nm,bio_index=di,n=nrow(dat),p=ncol(dat),replicate=rep,seed=seed,method=m,
        runtime_sec=elapsed,nrmse=met["nrmse"],pfc=met["pfc"],coefficient_abs_bias=bias,
        attempted=attempted,failed=nzchar(msg),failure=msg))
      write.csv(res,path,row.names=FALSE,na="")
    }
    message(nm," replicate ",rep,"/30")
  }
}
# Four metric wins (ties share a win) and mean dataset rank; NA metrics are omitted.
okres <- res[res$attempted & !res$failed, , drop=FALSE]
cell_groups <- split(okres, interaction(okres$dataset,okres$method,drop=TRUE))
mean_or_na <- function(x) if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)
agg <- do.call(rbind,lapply(cell_groups,function(d) data.frame(
  dataset=d$dataset[1],bio_index=d$bio_index[1],method=d$method[1],
  runtime_sec=mean_or_na(d$runtime_sec),nrmse=mean_or_na(d$nrmse),
  pfc=mean_or_na(d$pfc),coefficient_abs_bias=mean_or_na(d$coefficient_abs_bias))))
ms <- c("runtime_sec","nrmse","pfc","coefficient_abs_bias")
score <- lapply(split(agg,agg$dataset),function(d) {
  kk <- d[d$method=="missknn",,drop=FALSE]; ranks <- wins <- numeric()
  for(m in ms) {v<-d[[m]]; ok<-is.finite(v); if(!is.finite(kk[[m]][1])||!any(ok)) next
    ranks <- c(ranks,rank(v,ties.method="average",na.last="keep")[which(d$method=="missknn")])
    wins <- c(wins,as.numeric(kk[[m]][1]==min(v[ok]))) }
  data.frame(dataset=kk$dataset,bio_index=kk$bio_index,metric_wins=sum(wins),
    average_rank=if(length(ranks)) mean(ranks) else NA_real_)
})
score <- do.call(rbind,score); score <- score[order(-score$metric_wins,score$average_rank,score$bio_index),]
write.csv(agg,file.path(out,"biobench_dataset_summary.csv"),row.names=FALSE)
write.csv(score,file.path(out,"biobench_selection_ranking.csv"),row.names=FALSE)
write.csv(head(score,5),file.path(out,"biobench_selected_five.csv"),row.names=FALSE)
report <- c(
  "# Secondary BioBench regression benchmark", "",
  "Nine complete-case regression datasets; 30 replicates at 30% MCAR; four methods; one thread.",
  "Runtime, pooled NRMSE/PFC and downstream known-beta coefficient absolute bias are available in",
  "`biobench_all_results.csv`. The CSV records the seed and error text for every attempted cell.",
  "Failed cells are excluded from dataset means and rankings, and remain visible in the raw file.",
  "",
  paste0("Cells represented: ",nrow(res),"; attempted: ",sum(res$attempted),
         "; failed attempts: ",sum(res$attempted & res$failed),
         "; not attempted: ",sum(!res$attempted),"."),
  paste0("Results-selected five (metric wins, mean rank, BioBench order): ",
         paste(head(score$dataset,5),collapse=", "),"."),
  "",
  "See `biobench_selection_ranking.csv` for all nine datasets and the generated order, and",
  "`biobench_selected_five.csv` for the selected subset.",
  ""
)
writeLines(report,file.path(out,"biobench_report.md"))
