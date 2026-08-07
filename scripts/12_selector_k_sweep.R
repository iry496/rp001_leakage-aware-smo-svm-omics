#!/usr/bin/env Rscript
# =============================================================================
# scripts/12_selector_k_sweep.R
# -----------------------------------------------------------------------------
# Sensitivity analysis: vary the t-test/Welch top-K feature budget and compare
# the leaky vs guarded nested SMO/SVM pipelines (GSE25055 ONLY).
#
# Selector is held fixed (Welch two-sample t-test, top-K); only K varies over
# {25, 50, 100, 200}. For each K we compute, at the project seed (20260620):
#   * leaky AUROC / PR-AUC      (global FS top-K before CV, 5x5 repeated CV, cost 1)
#   * guarded AUROC / PR-AUC    (nested 5-outer x 5-inner; FS + cost tuning in folds)
#   * descriptive whole-workflow contrast dAUROC / dPR-AUC (naive - guarded)
#   * guarded feature stability: Nogueira index, mean Jaccard, stable-core and
#     unstable-tail counts (from the five outer-fold selected-feature sets)
#
# K = 100 is the ANCHOR and must reproduce the committed pilot:
#   naive AUROC ~0.7830, guarded ~0.7032, Nogueira ~0.5128, mean Jaccard ~0.3487,
#   stable core = 26, unstable tail = 105. If not, the run STOPS.
#
# Reuses the validated, optimized machinery from scripts 08/09 (fast Welch
# selector + selector validation guard + rds cache). NO external cohorts.
#
# Run (auto-runs main() via the sys.nframe guard):
#   Rscript scripts/12_selector_k_sweep.R smoke   # K=100 anchor only (fast check)
#   Rscript scripts/12_selector_k_sweep.R full     # K = 25, 50, 100, 200
# =============================================================================

suppressWarnings(suppressMessages({
  source("R/00_config.R")          # SEED (20260620)
  source("R/feature_selection.R")  # select_features() (original t.test selector)
  source("R/preprocessing.R")      # filter_near_zero_variance(), standardizer
  source("R/model_smo_svm.R")      # train_smo_svm(), predict_smo_svm()
  source("R/metrics.R")            # compute_binary_metrics()
  source("R/leakage_checks.R")     # assert_no_overlap()
}))

# ---- config -----------------------------------------------------------------
RUN_TAG <- "full"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1 && nzchar(.args[1])) RUN_TAG <- .args[1]
K_GRID  <- if (identical(RUN_TAG, "smoke")) c(100L) else c(25L, 50L, 100L, 200L)
USE_FAST_SELECTOR <- TRUE

ACCESSION <- "GSE25055"; LABEL_FIELD <- "pathologic_response_pcr_rd"; POSITIVE <- "pCR"
KERNEL <- "linear"
N_FOLDS <- 5; N_REPEATS <- 5; COST <- 1                       # leaky arm
OUTER_FOLDS <- 5; INNER_FOLDS <- 5; COST_GRID <- c(0.25, 1, 4)  # guarded arm
# anchor (K=100) committed references
REF <- list(leaky_auroc = 0.7830, guarded_auroc = 0.7032,
            nogueira = 0.5128, mean_jaccard = 0.3487, core = 26L, tail = 105L)

TAB_DIR <- file.path("tables", "sensitivity")
RES_DIR <- file.path("results", "sensitivity")
FIG_DIR <- file.path("figures", "sensitivity")
CACHE_X <- file.path("processed_data", "gse25055_perm_cache.rds")   # gitignored (shared w/ 08/09)
CACHE_Y <- file.path("processed_data", "gse25055_perm_labels.rds")
for (d in c(TAB_DIR, RES_DIR, FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- data loader (verbatim from scripts 08/09) ------------------------------
load_gse25055 <- function(accession = ACCESSION, label_field = LABEL_FIELD) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery required.")
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase required.")
  gset <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE); eset <- gset[[1]]
  expr <- Biobase::exprs(eset); pheno <- Biobase::pData(eset)
  label_col <- NULL
  for (cn in colnames(pheno)) {
    vals <- as.character(pheno[[cn]])
    if (any(grepl(label_field, vals, ignore.case = TRUE)) || grepl(label_field, cn, ignore.case = TRUE)) { label_col <- cn; break }
  }
  if (is.null(label_col)) stop("Label field not found.")
  raw <- as.character(pheno[[label_col]])
  raw <- trimws(sub(paste0("^.*", label_field, "\\s*[:=]?\\s*"), "", raw, ignore.case = TRUE))
  keep <- !(is.na(raw) | raw %in% c("NA","na","N/A","","NaN"))
  labels <- factor(raw[keep], levels = c("RD", "pCR"))
  x <- t(expr[, keep, drop = FALSE]); rownames(x) <- colnames(expr)[keep]
  ord <- order(rownames(x)); x <- x[ord, , drop = FALSE]; labels <- labels[ord]
  list(x = x, y = labels)
}
prep_expression <- function(x) filter_near_zero_variance(x, cutoff = 1e-8)$x
load_or_cache <- function() {
  if (file.exists(CACHE_X) && file.exists(CACHE_Y)) {
    message("[cache] reading processed matrix/labels from rds cache.")
    x <- readRDS(CACHE_X); y <- readRDS(CACHE_Y)
    ord <- order(rownames(x))
    return(list(x = x[ord, , drop = FALSE], y = y[ord]))
  }
  dat <- load_gse25055(); x <- prep_expression(dat$x); y <- dat$y
  dir.create("processed_data", recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, CACHE_X); saveRDS(y, CACHE_Y); list(x = x, y = y)
}

# ---- fast Welch selector (verbatim from 08/09) + dispatcher + validation -----
select_ttest_fast <- function(x_train, y_train, top_k) {
  x <- as.matrix(x_train); y <- as.factor(y_train); lv <- levels(y)
  g1 <- which(y == lv[1]); g2 <- which(y == lv[2]); n1 <- length(g1); n2 <- length(g2)
  m1 <- colMeans(x[g1,,drop=FALSE]); m2 <- colMeans(x[g2,,drop=FALSE])
  v1 <- (colSums(x[g1,,drop=FALSE]^2) - n1*m1^2)/(n1-1); v2 <- (colSums(x[g2,,drop=FALSE]^2) - n2*m2^2)/(n2-1)
  v1[v1<0] <- 0; v2[v2<0] <- 0
  s <- abs((m1-m2)/sqrt(v1/n1+v2/n2)); s[!is.finite(s)] <- 0; names(s) <- colnames(x)
  names(sort(s, decreasing = TRUE))[seq_len(min(top_k, length(s)))]
}
fs_select <- function(x_train, y_train, top_k) {
  if (isTRUE(USE_FAST_SELECTOR)) select_ttest_fast(x_train, y_train, top_k)
  else select_features(x_train, y_train, method = "t_test", top_k = top_k)$features
}
validate_selector <- function(x, y, top_k, n_check = 3) {
  cases <- c(list(identity = y), setNames(lapply(seq_len(n_check), function(b){set.seed(SEED+900000+b); y[sample.int(length(y))]}), paste0("perm", seq_len(n_check))))
  do.call(rbind, lapply(names(cases), function(nm){
    yy <- cases[[nm]]
    old <- select_features(x, yy, method="t_test", top_k=top_k)$features
    new <- select_ttest_fast(x, yy, top_k)
    data.frame(case=nm, top_k=top_k, overlap=length(intersect(old,new)), exact_order=identical(old,new), stringsAsFactors=FALSE)
  }))
}

# ---- leaky arm (logic from scripts 02/08; seed + top_k) ---------------------
run_leaky <- function(x, y, top_k, seed = SEED) {
  sel <- fs_select(x, y, top_k); x_sel <- x[, sel, drop=FALSE]; preds <- list()
  for (rep in seq_len(N_REPEATS)) {
    set.seed(seed + rep); folds <- caret::createFolds(y, k=N_FOLDS, returnTrain=FALSE)
    for (f in seq_along(folds)) {
      te <- folds[[f]]; tr <- setdiff(seq_along(y), te)
      std <- fit_standardizer(x_sel[tr,,drop=FALSE])
      model <- train_smo_svm(apply_standardizer(x_sel[tr,,drop=FALSE],std), y[tr], kernel=KERNEL, cost=COST, class_weights=TRUE, probability=TRUE)
      pr <- predict_smo_svm(model, apply_standardizer(x_sel[te,,drop=FALSE],std), positive_class=POSITIVE)
      preds[[length(preds)+1]] <- data.frame(repeat_id=rep, truth=as.character(y[te]), predicted=pr$predicted_class, prob_pos=pr$probability_positive, stringsAsFactors=FALSE)
    }
  }
  pd <- do.call(rbind, preds)
  agg <- do.call(rbind, lapply(split(pd, pd$repeat_id), function(d){
    m <- compute_binary_metrics(d$truth, d$predicted, d$prob_pos, POSITIVE)
    sens <- sum(d$truth==POSITIVE & d$predicted==POSITIVE)/sum(d$truth==POSITIVE)
    spec <- sum(d$truth!=POSITIVE & d$predicted!=POSITIVE)/sum(d$truth!=POSITIVE)
    data.frame(auroc=m$auroc, pr_auc=m$pr_auc, bal_acc=m$balanced_accuracy, mcc=m$mcc,
               sensitivity=sens, specificity=spec)}))
  c(auroc=mean(agg$auroc), pr_auc=mean(agg$pr_auc), bal_acc=mean(agg$bal_acc),
    mcc=mean(agg$mcc), sensitivity=mean(agg$sensitivity), specificity=mean(agg$specificity))
}

# ---- guarded nested arm (logic from scripts 03/09; seed + top_k) ------------
tune_inner <- function(x_tr, y_tr, top_k, seed = SEED) {
  set.seed(seed); folds <- caret::createFolds(y_tr, k=INNER_FOLDS, returnTrain=FALSE)
  score <- function(cost){ bal <- c()
    for (f in seq_along(folds)) {
      va <- folds[[f]]; tr <- setdiff(seq_along(y_tr), va); feats <- fs_select(x_tr[tr,,drop=FALSE], y_tr[tr], top_k)
      std <- fit_standardizer(x_tr[tr,feats,drop=FALSE])
      model <- train_smo_svm(apply_standardizer(x_tr[tr,feats,drop=FALSE],std), y_tr[tr], kernel=KERNEL, cost=cost, class_weights=TRUE, probability=TRUE)
      pr <- predict_smo_svm(model, apply_standardizer(x_tr[va,feats,drop=FALSE],std), positive_class=POSITIVE)
      sens <- sum(y_tr[va]==POSITIVE & pr$predicted_class==POSITIVE)/sum(y_tr[va]==POSITIVE)
      spec <- sum(y_tr[va]!=POSITIVE & pr$predicted_class!=POSITIVE)/sum(y_tr[va]!=POSITIVE)
      bal <- c(bal, mean(c(sens,spec), na.rm=TRUE)) }
    mean(bal, na.rm=TRUE) }
  COST_GRID[which.max(vapply(COST_GRID, score, numeric(1)))]
}
run_guarded <- function(x, y, top_k, seed = SEED) {
  set.seed(seed); folds <- caret::createFolds(y, k=OUTER_FOLDS, returnTrain=FALSE)
  preds <- list(); fold_feats <- vector("list", length(folds))
  for (f in seq_along(folds)) {
    te <- folds[[f]]; tr <- setdiff(seq_along(y), te)
    assert_no_overlap(rownames(x)[tr], rownames(x)[te])
    xt <- x[tr,,drop=FALSE]; yt <- y[tr]; bc <- tune_inner(xt, yt, top_k, seed)
    feats <- fs_select(xt, yt, top_k); fold_feats[[f]] <- feats
    std <- fit_standardizer(xt[,feats,drop=FALSE])
    model <- train_smo_svm(apply_standardizer(xt[,feats,drop=FALSE],std), yt, kernel=KERNEL, cost=bc, class_weights=TRUE, probability=TRUE)
    pr <- predict_smo_svm(model, apply_standardizer(x[te,feats,drop=FALSE],std), positive_class=POSITIVE)
    preds[[f]] <- data.frame(truth=as.character(y[te]), predicted=pr$predicted_class, prob_pos=pr$probability_positive, stringsAsFactors=FALSE)
  }
  pd <- do.call(rbind, preds); m <- compute_binary_metrics(pd$truth, pd$predicted, pd$prob_pos, POSITIVE)
  sens <- sum(pd$truth==POSITIVE & pd$predicted==POSITIVE)/sum(pd$truth==POSITIVE)
  spec <- sum(pd$truth!=POSITIVE & pd$predicted!=POSITIVE)/sum(pd$truth!=POSITIVE)
  list(metrics=c(auroc=m$auroc, pr_auc=m$pr_auc, bal_acc=m$balanced_accuracy, mcc=m$mcc,
                 sensitivity=sens, specificity=spec), fold_feats=fold_feats)
}
stability_from_feats <- function(fold_feats, P) {
  M <- length(fold_feats); k_each <- vapply(fold_feats, length, integer(1)); kbar <- mean(k_each)
  freq <- as.integer(table(unlist(fold_feats))); phat <- freq/M
  nog <- 1 - (sum((M/(M-1))*phat*(1-phat))/P) / ((kbar/P)*(1-kbar/P))
  pair <- c(); for (i in 1:(M-1)) for (j in (i+1):M){a<-fold_feats[[i]];b<-fold_feats[[j]];pair<-c(pair,length(intersect(a,b))/length(union(a,b)))}
  c(nogueira=nog, mean_jaccard=mean(pair), total_unique=length(freq), stable_core=sum(freq==M), unstable_tail=sum(freq==1))
}

# ---- driver -----------------------------------------------------------------
main <- function() {
  t0 <- Sys.time()
  message(sprintf("[cfg] RUN_TAG=%s K_GRID=%s FAST=%s", RUN_TAG, paste(K_GRID, collapse=","), USE_FAST_SELECTOR))
  dd <- load_or_cache(); x <- dd$x; y <- dd$y; P <- ncol(x)
  message(sprintf("[ksweep] N=%d; P=%d; class %s", length(y), P, paste(names(table(y)), table(y), sep="=", collapse=", ")))
  if (isTRUE(USE_FAST_SELECTOR)) { message("[validate] selector check at K=100 ..."); print(validate_selector(x, y, 100L, 3)) }

  rows <- list()
  for (K in K_GRID) {
    tk <- Sys.time(); lk <- run_leaky(x, y, K); gd <- run_guarded(x, y, K)
    st <- stability_from_feats(gd$fold_feats, P)
    rows[[length(rows)+1]] <- data.frame(
      top_k=K,
      leaky_auroc=round(lk["auroc"],4), guarded_auroc=round(gd$metrics["auroc"],4),
      gap_auroc=round(lk["auroc"]-gd$metrics["auroc"],4),
      leaky_pr_auc=round(lk["pr_auc"],4), guarded_pr_auc=round(gd$metrics["pr_auc"],4),
      gap_pr_auc=round(lk["pr_auc"]-gd$metrics["pr_auc"],4),
      leaky_bal_acc=round(lk["bal_acc"],4), guarded_bal_acc=round(gd$metrics["bal_acc"],4),
      leaky_mcc=round(lk["mcc"],4), guarded_mcc=round(gd$metrics["mcc"],4),
      leaky_sensitivity=round(lk["sensitivity"],4), guarded_sensitivity=round(gd$metrics["sensitivity"],4),
      leaky_specificity=round(lk["specificity"],4), guarded_specificity=round(gd$metrics["specificity"],4),
      nogueira=round(st["nogueira"],4), mean_jaccard=round(st["mean_jaccard"],4),
      total_unique=st["total_unique"], stable_core=st["stable_core"], unstable_tail=st["unstable_tail"],
      runtime_secs=round(as.numeric(difftime(Sys.time(),tk,units="secs")),1), stringsAsFactors=FALSE)
    message(sprintf("[K=%d] leaky=%.4f guarded=%.4f gap=%.4f | Nogueira=%.4f Jaccard=%.4f core=%d tail=%d (%.1fs)",
                    K, lk["auroc"], gd$metrics["auroc"], lk["auroc"]-gd$metrics["auroc"], st["nogueira"], st["mean_jaccard"], st["stable_core"], st["unstable_tail"],
                    as.numeric(difftime(Sys.time(),tk,units="secs"))))
  }
  res <- do.call(rbind, rows); rownames(res) <- NULL
  write.csv(res, file.path(TAB_DIR, "k_sweep_summary.csv"), row.names = FALSE)

  # anchor reproduction check at K=100
  a <- res[res$top_k==100L, ]
  ok <- TRUE
  if (nrow(a)==1) {
    chk <- c(abs(a$leaky_auroc-REF$leaky_auroc)<0.01, abs(a$guarded_auroc-REF$guarded_auroc)<0.01,
             abs(a$nogueira-REF$nogueira)<0.01, abs(a$mean_jaccard-REF$mean_jaccard)<0.01,
             a$stable_core==REF$core, a$unstable_tail==REF$tail)
    ok <- all(chk)
    message(sprintf("[anchor K=100] leaky %.4f/%.4f guarded %.4f/%.4f Nogueira %.4f/%.4f Jaccard %.4f/%.4f core %d/%d tail %d/%d -> %s",
                    a$leaky_auroc,REF$leaky_auroc,a$guarded_auroc,REF$guarded_auroc,a$nogueira,REF$nogueira,
                    a$mean_jaccard,REF$mean_jaccard,a$stable_core,REF$core,a$unstable_tail,REF$tail, if(ok)"PASS" else "REVIEW"))
    if (!ok) stop("[ksweep] anchor K=100 did not reproduce committed pilot; STOP and debug before trusting the sweep.")
  }

  # figures (only for full sweep with >1 K)
  if (length(K_GRID) > 1) {
    Kx <- res$top_k
    for (dev_fun in list(function() png(file.path(FIG_DIR,"k_sweep_leakage_gap.png"), width=1100, height=500),
                         function() pdf(file.path(FIG_DIR,"k_sweep_leakage_gap.pdf"), width=11, height=5))) {
      dev_fun(); par(mfrow=c(1,2))
      plot(Kx, res$leaky_auroc, type="b", pch=19, col="#C0392B", ylim=range(c(res$leaky_auroc,res$guarded_auroc)),
           xlab="top-K features", ylab="AUROC", main="Leaky vs guarded AUROC by K", xaxt="n"); axis(1, at=Kx)
      lines(Kx, res$guarded_auroc, type="b", pch=17, col="#1E8449"); legend("bottomright", pch=c(19,17), col=c("#C0392B","#1E8449"), legend=c("leaky","guarded"), bty="n")
      plot(Kx, res$gap_auroc, type="b", pch=19, col="#34495E", ylim=range(c(0,res$gap_auroc,res$gap_pr_auc)),
           xlab="top-K features", ylab="whole-workflow contrast", main="Naive - guarded contrast by K", xaxt="n"); axis(1, at=Kx)
      lines(Kx, res$gap_pr_auc, type="b", pch=15, col="#7D6608"); abline(h=0, lty=3, col="#888888")
      legend("topright", pch=c(19,15), col=c("#34495E","#7D6608"), legend=c("dAUROC","dPR-AUC"), bty="n"); dev.off()
    }
    for (dev_fun in list(function() png(file.path(FIG_DIR,"k_sweep_stability.png"), width=1100, height=500),
                         function() pdf(file.path(FIG_DIR,"k_sweep_stability.pdf"), width=11, height=5))) {
      dev_fun(); par(mfrow=c(1,2))
      plot(Kx, res$nogueira, type="b", pch=19, col="#6C3483", ylim=range(c(res$nogueira,res$mean_jaccard)),
           xlab="top-K features", ylab="stability index", main="Stability by K", xaxt="n"); axis(1, at=Kx)
      lines(Kx, res$mean_jaccard, type="b", pch=17, col="#1F618D"); legend("topright", pch=c(19,17), col=c("#6C3483","#1F618D"), legend=c("Nogueira","mean Jaccard"), bty="n")
      plot(Kx, res$stable_core, type="b", pch=19, col="#1E8449", ylim=range(c(res$stable_core,res$unstable_tail)),
           xlab="top-K features", ylab="probe count", main="Stable core vs unstable tail by K", xaxt="n"); axis(1, at=Kx)
      lines(Kx, res$unstable_tail, type="b", pch=17, col="#C0392B"); legend("topleft", pch=c(19,17), col=c("#1E8449","#C0392B"), legend=c("stable core","unstable tail"), bty="n"); dev.off()
    }
  }

  notes <- c(
    sprintf("# Selector K-sweep (%s run, GSE25055 only)", RUN_TAG), "",
    sprintf("- Selector fixed (Welch t-test top-K); K in {%s}. Seed %d.", paste(K_GRID, collapse=", "), SEED),
    "- Leaky: global FS top-K before 5x5 repeated CV (cost 1). Guarded: nested 5-outer x 5-inner; FS + cost tuning in training folds.",
    "- K=100 is the anchor and reproduces the submission-aligned run (naive ~0.7830, guarded ~0.7032, Nogueira ~0.5128).",
    "", "## Results by K")
  for (i in seq_len(nrow(res))) notes <- c(notes, sprintf(
    "- K=%d: AUROC leaky %.4f / guarded %.4f (gap %.4f); PR-AUC leaky %.4f / guarded %.4f (gap %.4f); bal.acc %.4f / %.4f; MCC %.4f / %.4f; sens %.4f / %.4f; spec %.4f / %.4f; Nogueira %.4f, mean Jaccard %.4f, stable core %d, unstable tail %d.",
    res$top_k[i], res$leaky_auroc[i], res$guarded_auroc[i], res$gap_auroc[i],
    res$leaky_pr_auc[i], res$guarded_pr_auc[i], res$gap_pr_auc[i],
    res$leaky_bal_acc[i], res$guarded_bal_acc[i], res$leaky_mcc[i], res$guarded_mcc[i],
    res$leaky_sensitivity[i], res$guarded_sensitivity[i], res$leaky_specificity[i], res$guarded_specificity[i],
    res$nogueira[i], res$mean_jaccard[i], res$stable_core[i], res$unstable_tail[i]))
  notes <- c(notes, "", sprintf("- Total runtime: %.1f min.", as.numeric(difftime(Sys.time(),t0,units="mins"))),
             "", "_Sensitivity analysis only; within-cohort, diagnostic; not biomarker discovery._")
  writeLines(notes, file.path(RES_DIR, sprintf("k_sweep_notes%s.md", if (identical(RUN_TAG,"smoke")) "_smoke" else "")))
  message(sprintf("[ksweep] DONE (%s). Wrote k_sweep_summary.csv%s.", RUN_TAG, if (length(K_GRID)>1) " + figures" else " (anchor smoke; no figures)"))
  invisible(res)
}

if (sys.nframe() == 0) main()
