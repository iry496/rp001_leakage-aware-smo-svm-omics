#!/usr/bin/env Rscript
# =============================================================================
# scripts/14_external_validation_gse41998.R
# -----------------------------------------------------------------------------
# Cross-platform external validation of the FROZEN GSE25055 model on GSE41998,
# per docs/gse41998_external_validation_design.md (approved).
#
# LEAKAGE SAFE — GSE41998 is used ONLY for final evaluation:
#   * Feature set, SVM cost, scaling parameters, and the 0.5 P(pCR) threshold are
#     FROZEN on GSE25055 discovery; nothing is tuned/selected on GSE41998.
#   * Exact probe-ID intersection (GPL96 model features INTERSECT GPL571 probes);
#     no gene-symbol collapse here; no joint normalization; no ComBat.
#   * Predeclared label rule: include unambiguous pCR(Yes)/RD(No) only; exclude
#     "0", missing, and unclear labels.
#   * Scaling: PRIMARY = discovery-derived (frozen); SENSITIVITY = within-cohort
#     label-blind z-score.
#
# RUNTIME GO/NO-GO (from the memo) — the script STOPS (report-only) if not met:
#   * usable labeled n >= MIN_N with BOTH classes and non-extreme prevalence;
#   * probe intersection recovers >= MIN_COVERAGE of the frozen feature set.
#
# Does NOT write raw expression. Run: Rscript scripts/14_external_validation_gse41998.R
# =============================================================================

suppressWarnings(suppressMessages({
  source("R/00_config.R"); source("R/feature_selection.R"); source("R/preprocessing.R")
  source("R/model_smo_svm.R"); source("R/metrics.R"); source("R/leakage_checks.R")
}))

# ---- config (all frozen / predeclared) --------------------------------------
DISC <- "GSE25055"; EXT <- "GSE41998"; LABEL_FIELD <- "pathologic_response_pcr_rd"
POSITIVE <- "pCR"; TOP_K <- 100L; COST_GRID <- c(0.25, 1, 4); INNER_FOLDS <- 5; KERNEL <- "linear"
THRESHOLD <- 0.5                 # on P(pCR); no post-hoc tuning
MIN_N <- 80L; MIN_COVERAGE <- 0.70   # go/no-go thresholds (memo)
# GSE41998 label tokens (predeclared); positive = pCR/Yes, negative = RD/No; everything else excluded
POS_TOK <- c("pcr","yes","y","1","complete response","pathologic complete response")
NEG_TOK <- c("rd","no","n","residual","residual disease","non-pcr","not pcr")
DISC_CACHE_X <- file.path("processed_data","gse25055_perm_cache.rds")
DISC_CACHE_Y <- file.path("processed_data","gse25055_perm_labels.rds")
TAB <- "tables/external_validation_gse41998"; RES <- "results/external_validation_gse41998"; FIG <- "figures/external_validation_gse41998"
for (d in c(TAB,RES,FIG)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

stop_report_only <- function(reason, extra=list()) {
  notes <- c("# GSE41998 external validation — NO-GO (report only)", "",
             paste0("Status: NOT RUN. Reason: ", reason), "",
             "Per the approved design memo, the frozen-model projection was not produced because a",
             "go/no-go criterion was not met. No metrics are reported. GSE41998 remains future work.")
  for (k in names(extra)) notes <- c(notes, sprintf("- %s: %s", k, extra[[k]]))
  writeLines(notes, file.path(RES, "gse41998_notes.md"))
  write.csv(data.frame(status="no_go", reason=reason, stringsAsFactors=FALSE),
            file.path(TAB, "gse41998_summary.csv"), row.names=FALSE)
  message("[gse41998] NO-GO: ", reason, " (wrote report-only notes/summary)."); quit(save="no", status=0)
}

# ---- discovery (cache shared w/ 08/09) --------------------------------------
load_discovery <- function() {
  if (file.exists(DISC_CACHE_X) && file.exists(DISC_CACHE_Y)) {
    message("[cache] reading GSE25055 discovery matrix/labels from rds cache.")
    x <- readRDS(DISC_CACHE_X); y <- readRDS(DISC_CACHE_Y)
    ord <- order(rownames(x))
    return(list(x=x[ord, , drop=FALSE], y=y[ord]))
  }
  if (!requireNamespace("GEOquery", quietly=TRUE)) stop("GEOquery required.")
  g <- GEOquery::getGEO(DISC, GSEMatrix=TRUE, getGPL=FALSE)[[1]]
  expr <- Biobase::exprs(g); ph <- Biobase::pData(g); lc <- NULL
  for (cn in colnames(ph)) if (any(grepl(LABEL_FIELD, as.character(ph[[cn]]), ignore.case=TRUE)) || grepl(LABEL_FIELD, cn, ignore.case=TRUE)) {lc<-cn;break}
  raw <- trimws(sub(paste0("^.*",LABEL_FIELD,"\\s*[:=]?\\s*"),"",as.character(ph[[lc]]),ignore.case=TRUE))
  keep <- !(is.na(raw)|raw %in% c("NA","na","N/A","","NaN")); y <- factor(raw[keep], levels=c("RD","pCR"))
  x <- t(expr[,keep,drop=FALSE]); rownames(x) <- colnames(expr)[keep]
  ord <- order(rownames(x)); x <- x[ord, , drop=FALSE]; y <- y[ord]
  x <- filter_near_zero_variance(x, cutoff=1e-8)$x; list(x=x, y=y)
}

# ---- GSE41998 loader + predeclared label parsing ---------------------------
map_label <- function(v) {
  s <- tolower(trimws(as.character(v)))
  if (s %in% c("","na","n/a","nan","unknown","unclear","0","missing")) return(NA_character_)
  if (any(vapply(POS_TOK, function(t) grepl(paste0("(^|[^a-z])",t,"([^a-z]|$)"), s), logical(1)))) return("pCR")
  if (any(vapply(NEG_TOK, function(t) grepl(paste0("(^|[^a-z])",t,"([^a-z]|$)"), s), logical(1)))) return("RD")
  NA_character_
}
load_ext <- function() {
  if (!requireNamespace("GEOquery", quietly=TRUE)) stop("GEOquery required.")
  g <- GEOquery::getGEO(EXT, GSEMatrix=TRUE, getGPL=FALSE)
  g <- g[[1]]; expr <- Biobase::exprs(g); ph <- Biobase::pData(g)
  # PREDECLARED label resolution (per design memo): GSE41998 documents the binary endpoint in the
  # "pcr" characteristic, with values "pcr: Yes" / "pcr: No" / "pcr: 0". Select that EXACT field
  # (value prefix "pcr:"; this excludes "pcrrcb1:"), strip the field-name prefix, then map strictly
  # Yes -> pCR and No -> RD. The "0" code and missing/blank are not-evaluable and are EXCLUDED.
  # No GSE41998 performance is consulted; this reads the documented endpoint only
  # (expected 69 pCR / 184 RD = 253 evaluable). Earlier auto-detection failed because it applied the
  # token mapper to the un-stripped "pcr: ..." string (the field name "pcr" matched POS_TOK), so every
  # value mapped to pCR; stripping the prefix and mapping Yes/No strictly fixes that.
  is_pcr_col <- function(cn) {
    v <- tolower(trimws(as.character(ph[[cn]])))
    mean(grepl("^pcr\\s*:", v)) > 0.5
  }
  cand <- colnames(ph)[vapply(colnames(ph), is_pcr_col, logical(1))]
  if (!length(cand))
    stop_report_only("GSE41998 'pcr' label field (values 'pcr: Yes/No/0') not found in series-matrix metadata")
  lc <- cand[1]
  bare <- tolower(trimws(sub("^\\s*pcr\\s*:\\s*", "", as.character(ph[[lc]]), ignore.case=TRUE)))
  y_all <- ifelse(bare %in% c("yes","y"), "pCR",
           ifelse(bare %in% c("no","n"),  "RD", NA_character_))   # "0"/blank/NA/other -> excluded
  message(sprintf("[ext] label column '%s'; Yes->pCR / No->RD; excluded (0/missing/other)=%d of %d",
                  lc, sum(is.na(y_all)), nrow(ph)))
  y_all <- factor(y_all, levels=c("RD","pCR"))
  list(expr=expr, y=y_all, label_col=lc)
}

main <- function() {
  t0 <- Sys.time()
  dd <- load_discovery(); xd <- dd$x; yd <- dd$y
  message(sprintf("[disc] N=%d P=%d class %s", length(yd), ncol(xd), paste(names(table(yd)),table(yd),sep="=",collapse=", ")))

  # ---- FROZEN model on discovery (feature set, cost, scaling, threshold) ----
  feats <- select_features(xd, yd, method="t_test", top_k=TOP_K)$features   # top-100 on full discovery
  # cost by guarded discovery-only CV (no GSE41998)
  set.seed(SEED); folds <- caret::createFolds(yd, k=INNER_FOLDS, returnTrain=FALSE)
  score <- function(cost){ b<-c(); for(f in seq_along(folds)){va<-folds[[f]];tr<-setdiff(seq_along(yd),va)
      ff<-select_features(xd[tr,,drop=FALSE],yd[tr],method="t_test",top_k=TOP_K)$features
      std<-fit_standardizer(xd[tr,ff,drop=FALSE]); m<-train_smo_svm(apply_standardizer(xd[tr,ff,drop=FALSE],std),yd[tr],kernel=KERNEL,cost=cost,class_weights=TRUE,probability=TRUE)
      pr<-predict_smo_svm(m,apply_standardizer(xd[va,ff,drop=FALSE],std),positive_class=POSITIVE)
      se<-sum(yd[va]==POSITIVE&pr$predicted_class==POSITIVE)/sum(yd[va]==POSITIVE); sp<-sum(yd[va]!=POSITIVE&pr$predicted_class!=POSITIVE)/sum(yd[va]!=POSITIVE)
      b<-c(b,mean(c(se,sp),na.rm=TRUE))}; mean(b,na.rm=TRUE) }
  best_cost <- COST_GRID[which.max(vapply(COST_GRID, score, numeric(1)))]
  message(sprintf("[frozen] top-%d features; cost=%s (guarded discovery CV); threshold=%.2f", TOP_K, best_cost, THRESHOLD))

  # ---- load GSE41998 + predeclared labels ----------------------------------
  ee <- load_ext(); expr_e <- ee$expr; y_e_all <- ee$y
  keep <- !is.na(y_e_all)
  y_e <- y_e_all[keep]; expr_e <- expr_e[, keep, drop=FALSE]
  ncl <- table(y_e); prev <- as.numeric(ncl["pCR"]) / sum(ncl)
  message(sprintf("[ext] usable n=%d; class %s; pCR prevalence %.3f", length(y_e), paste(names(ncl),ncl,sep="=",collapse=", "), prev))

  # ---- GO/NO-GO: labels ----------------------------------------------------
  if (length(y_e) < MIN_N) stop_report_only(sprintf("usable labeled n=%d < MIN_N=%d", length(y_e), MIN_N))
  if (any(ncl < 5) || prev < 0.05 || prev > 0.60)
    stop_report_only(sprintf("class balance unsuitable (pCR=%s, RD=%s, prevalence=%.3f)", ncl["pCR"], ncl["RD"], prev))

  # ---- exact probe intersection (GPL96 model features INTERSECT GPL571) -----
  ext_probes <- rownames(expr_e)
  final_feats <- intersect(feats, ext_probes)
  coverage <- length(final_feats) / length(feats)
  message(sprintf("[map] frozen features recovered in GSE41998: %d/%d (coverage %.3f)", length(final_feats), length(feats), coverage))
  if (coverage < MIN_COVERAGE)
    stop_report_only(sprintf("probe-intersection coverage %.3f < MIN_COVERAGE=%.2f", coverage, MIN_COVERAGE),
                     list(recovered=length(final_feats), model_features=length(feats)))
  write.csv(data.frame(model_feature=feats, in_gse41998=feats %in% ext_probes),
            file.path(RES, "gse41998_probe_mapping.csv"), row.names=FALSE)

  # ---- refit frozen pipeline on the intersected feature set (discovery only;
  #      platform mapping fixed BEFORE using any GSE41998 label) --------------
  std_disc <- fit_standardizer(xd[, final_feats, drop=FALSE])              # discovery-derived scaling
  model    <- train_smo_svm(apply_standardizer(xd[, final_feats, drop=FALSE], std_disc), yd,
                            kernel=KERNEL, cost=best_cost, class_weights=TRUE, probability=TRUE)
  x_e <- t(expr_e[final_feats, , drop=FALSE]); rownames(x_e) <- colnames(expr_e)

  metr <- function(prob) {
    pred <- ifelse(prob >= THRESHOLD, POSITIVE, "RD")
    m <- compute_binary_metrics(as.character(y_e), pred, prob, POSITIVE)
    se <- sum(y_e==POSITIVE & pred==POSITIVE)/sum(y_e==POSITIVE); sp <- sum(y_e!=POSITIVE & pred!=POSITIVE)/sum(y_e!=POSITIVE)
    list(metrics=c(auroc=m$auroc, pr_auc=m$pr_auc, balanced_accuracy=m$balanced_accuracy, mcc=m$mcc, sensitivity=se, specificity=sp), pred=pred, prob=prob)
  }
  # PRIMARY: discovery-derived scaling
  pr_p <- predict_smo_svm(model, apply_standardizer(x_e[, final_feats, drop=FALSE], std_disc), positive_class=POSITIVE)
  res_p <- metr(pr_p$probability_positive)
  # SENSITIVITY: within-cohort label-blind z-score
  std_ext <- fit_standardizer(x_e[, final_feats, drop=FALSE])              # uses GSE41998 features only, NO labels
  pr_s <- predict_smo_svm(model, apply_standardizer(x_e[, final_feats, drop=FALSE], std_ext), positive_class=POSITIVE)
  res_s <- metr(pr_s$probability_positive)

  # ---- outputs -------------------------------------------------------------
  preds <- data.frame(sample_id=rownames(x_e), truth=as.character(y_e),
                      prob_pos_primary=round(res_p$prob,4), pred_primary=res_p$pred,
                      prob_pos_sensitivity=round(res_s$prob,4), pred_sensitivity=res_s$pred, stringsAsFactors=FALSE)
  write.csv(preds, file.path(RES, "gse41998_predictions.csv"), row.names=FALSE)
  mk <- function(tag, r) data.frame(scaling=tag, t(round(r$metrics,4)), check.names=FALSE, stringsAsFactors=FALSE)
  metrics_df <- rbind(mk("discovery_derived_primary", res_p), mk("within_cohort_zscore_sensitivity", res_s))
  write.csv(metrics_df, file.path(RES, "gse41998_metrics.csv"), row.names=FALSE)
  summ <- data.frame(status="go", external_cohort=EXT, usable_n=length(y_e), pCR=as.integer(ncl["pCR"]), RD=as.integer(ncl["RD"]),
                     prevalence=round(prev,3), feature_coverage=round(coverage,3), frozen_cost=best_cost, threshold=THRESHOLD,
                     auroc_primary=res_p$metrics["auroc"], pr_auc_primary=res_p$metrics["pr_auc"],
                     auroc_sensitivity=res_s$metrics["auroc"], pr_auc_sensitivity=res_s$metrics["pr_auc"], stringsAsFactors=FALSE)
  rownames(summ) <- NULL
  write.csv(summ, file.path(TAB, "gse41998_summary.csv"), row.names=FALSE)

  # ---- transportability figure (discovery guarded CV vs GSE25065 vs GSE41998) ----
  disc_metrics <- read.csv("results/pilot_gse25055/nested_smo_svm_metrics.csv", stringsAsFactors=FALSE)
  g65_metrics <- read.csv("results/external_validation_gse25065/gse25065_external_metrics.csv", stringsAsFactors=FALSE)
  disc_auroc <- disc_metrics$auroc[1]; disc_pr <- disc_metrics$pr_auc[1]
  g65_auroc <- g65_metrics$auroc[1]; g65_pr <- g65_metrics$pr_auc[1]
  M <- rbind(AUROC = c(disc_auroc, g65_auroc, res_p$metrics["auroc"]),
             `PR-AUC` = c(disc_pr, g65_pr, res_p$metrics["pr_auc"]))
  for (dev_fun in list(function() png(file.path(FIG,"gse41998_transportability.png"), width=900, height=520),
                       function() pdf(file.path(FIG,"gse41998_transportability.pdf"), width=9, height=5.2))) {
    dev_fun(); par(mar=c(4,4,3,1))
    bp <- barplot(M, beside=TRUE, col=c("#1E8449","#7FB3A6"), ylim=c(0,0.9),
                  names.arg=c("discovery (guarded CV)","GSE25065 (same-family)","GSE41998 (cross-platform)"),
                  legend.text=c("AUROC","PR-AUC"), args.legend=list(x="topright",bty="n"),
                  main="Transportability of the frozen model (primary, discovery scaling)")
    text(bp, c(M)+0.03, sprintf("%.2f", c(M)), cex=0.7); dev.off()
  }

  notes <- c(sprintf("# GSE41998 external validation (GO) — frozen GSE25055 model"), "",
    sprintf("- Status: GO. Usable n=%d (pCR=%s, RD=%s; prevalence %.3f). Label column: '%s'.", length(y_e), ncl["pCR"], ncl["RD"], prev, ee$label_col),
    sprintf("- Frozen feature set: top-%d (Welch t-test on full GSE25055); recovered in GSE41998: %d/%d (coverage %.3f).", TOP_K, length(final_feats), length(feats), coverage),
    sprintf("- Frozen cost=%s (guarded discovery CV); threshold=%.2f on P(pCR); no GSE41998 tuning/selection.", best_cost, THRESHOLD),
    "- Exact probe-ID intersection; no gene-symbol collapse; no joint normalization; no ComBat.",
    "", "## Results (frozen-model projection)",
    sprintf("- PRIMARY (discovery-derived scaling): AUROC %.4f, PR-AUC %.4f, balanced acc %.4f, MCC %.4f, sens %.4f, spec %.4f.",
            res_p$metrics["auroc"],res_p$metrics["pr_auc"],res_p$metrics["balanced_accuracy"],res_p$metrics["mcc"],res_p$metrics["sensitivity"],res_p$metrics["specificity"]),
    sprintf("- SENSITIVITY (within-cohort z-score): AUROC %.4f, PR-AUC %.4f.", res_s$metrics["auroc"], res_s$metrics["pr_auc"]),
    "", "## Interpretation",
    "Cross-platform transportability of a frozen model; a generalization limit, not a leakage effect.",
    "Within-cohort/diagnostic; not biomarker discovery. GSE41998 labels used only for final evaluation.",
    "", sprintf("_Runtime %.1f min. No raw expression written._", as.numeric(difftime(Sys.time(),t0,units="mins"))))
  writeLines(notes, file.path(RES, "gse41998_notes.md"))
  message(sprintf("[gse41998] GO. PRIMARY AUROC %.4f / PR-AUC %.4f ; SENSITIVITY AUROC %.4f. Outputs written.",
                  res_p$metrics["auroc"], res_p$metrics["pr_auc"], res_s$metrics["auroc"]))
}
if (sys.nframe() == 0) main()
