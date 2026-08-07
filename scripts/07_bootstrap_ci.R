#!/usr/bin/env Rscript
# =============================================================================
# scripts/07_bootstrap_ci.R
# -----------------------------------------------------------------------------
# Uncertainty quantification for the leakage-aware SMO/SVM audit.
#
# This script performs NO new modeling. It RESAMPLES the already-committed
# out-of-fold predictions to attach uncertainty to the headline metrics:
#   * Stratified percentile bootstrap 95% CIs for AUROC, PR-AUC, balanced
#     accuracy, MCC, sensitivity, and specificity, for the leaky baseline and
#     guarded nested pipeline (GSE25055) and the external cohort (GSE25065).
#   * Paired stratified bootstrap CI + two-sided bootstrap p-value for the
#     leakage gap (Delta AUROC and Delta PR-AUC, leaky - guarded).
#   * DeLong test (per repeat) for the AUROC difference, leaky vs guarded.
#
# Inputs (read-only, already committed):
#   results/pilot_gse25055/leaky_baseline_predictions.csv   (5 repeats x 306)
#   results/pilot_gse25055/nested_smo_svm_predictions.csv   (306)
#   results/external_validation_gse25065/gse25065_external_predictions.csv (182)
#
# Outputs:
#   tables/uncertainty/bootstrap_ci.csv
#   tables/uncertainty/delta_auroc_prauc_ci.csv
#   results/uncertainty/delong_tests.csv
#
# Notes:
#   * The leaky baseline AUROC/PR-AUC point estimates are the MEAN across the
#     5 repeated-CV runs (matching the manuscript), so the bootstrap recomputes
#     that mean within each resample. Operating-point metrics pool the repeats.
#   * Point estimates and DeLong p-values are deterministic. Bootstrap CIs are
#     reproducible within R given set.seed() below; a Python cross-check
#     (scripts/checks/07_bootstrap_ci_crosscheck.py) generated the committed
#     CSVs and agrees with this script on point estimates and DeLong results.
# =============================================================================

suppressPackageStartupMessages({
  library(pROC)
  library(PRROC)
})

set.seed(20260620)
B   <- 2000
POS <- "pCR"
NEG <- "RD"

dir.create("tables/uncertainty",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/uncertainty", recursive = TRUE, showWarnings = FALSE)

# ---- metric helpers ---------------------------------------------------------
# Reuse the project metric definitions where available.
if (file.exists("R/metrics.R")) source("R/metrics.R")

auroc <- function(truth, prob) {
  as.numeric(pROC::auc(pROC::roc(response = truth, predictor = prob,
                                 levels = c(NEG, POS), direction = "<",
                                 quiet = TRUE)))
}
pr_auc <- function(truth, prob) {
  tb <- as.integer(truth == POS)
  PRROC::pr.curve(scores.class0 = prob[tb == 1],
                  scores.class1 = prob[tb == 0], curve = FALSE)$auc.integral
}
op_metrics <- function(truth, pred) {
  sens <- sum(truth == POS & pred == POS) / sum(truth == POS)
  spec <- sum(truth != POS & pred != POS) / sum(truth != POS)
  tp <- as.numeric(sum(truth == POS & pred == POS)); tn <- as.numeric(sum(truth != POS & pred != POS))
  fp <- as.numeric(sum(truth != POS & pred == POS)); fn <- as.numeric(sum(truth == POS & pred != POS))
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (den == 0) NA_real_ else (tp * tn - fp * fn) / den
  c(balanced_accuracy = mean(c(sens, spec)), mcc = mcc,
    sensitivity = sens, specificity = spec)
}
METRICS <- c("auroc", "pr_auc", "balanced_accuracy", "mcc",
             "sensitivity", "specificity")

# ---- load predictions -------------------------------------------------------
leaky  <- read.csv("results/pilot_gse25055/leaky_baseline_predictions.csv",
                   stringsAsFactors = FALSE)
nested <- read.csv("results/pilot_gse25055/nested_smo_svm_predictions.csv",
                   stringsAsFactors = FALSE)
ext    <- read.csv("results/external_validation_gse25065/gse25065_external_predictions.csv",
                   stringsAsFactors = FALSE)

reps <- sort(unique(leaky$repeat_id))
ids  <- nested$sample_id                                   # 306 internal samples
truth_int <- setNames(nested$truth, nested$sample_id)[ids]
# per-repeat lookup tables for the leaky baseline
leaky_prob <- lapply(reps, function(r) {
  s <- leaky[leaky$repeat_id == r, ]; setNames(s$prob_pos, s$sample_id)[ids]
})
nested_prob <- setNames(nested$prob_pos, nested$sample_id)[ids]
nested_pred <- setNames(nested$predicted, nested$sample_id)[ids]

# ---- point estimates --------------------------------------------------------
leaky_point <- c(
  auroc  = mean(vapply(leaky_prob, function(p) auroc(truth_int, p), numeric(1))),
  pr_auc = mean(vapply(leaky_prob, function(p) pr_auc(truth_int, p), numeric(1))),
  op_metrics(leaky$truth, leaky$predicted)                # pooled across repeats
)
nested_point <- c(auroc = auroc(truth_int, nested_prob),
                  pr_auc = pr_auc(truth_int, nested_prob),
                  op_metrics(truth_int, nested_pred))
ext_point <- c(auroc = auroc(ext$truth, ext$prob_pos),
               pr_auc = pr_auc(ext$truth, ext$prob_pos),
               op_metrics(ext$truth, ext$predicted))

# ---- stratified resampling helpers -----------------------------------------
pos_ids <- ids[truth_int == POS]; neg_ids <- ids[truth_int == NEG]
resample_int <- function() c(sample(pos_ids, replace = TRUE),
                             sample(neg_ids, replace = TRUE))
ext_pos <- which(ext$truth == POS); ext_neg <- which(ext$truth == NEG)
resample_ext <- function() c(sample(ext_pos, replace = TRUE),
                             sample(ext_neg, replace = TRUE))

pct_ci <- function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

# ---- bootstrap each cohort --------------------------------------------------
boot_cohort <- function(kind) {
  acc <- matrix(NA_real_, nrow = B, ncol = length(METRICS),
                dimnames = list(NULL, METRICS))
  for (b in seq_len(B)) {
    if (kind == "external") {
      idx <- resample_ext(); tr <- ext$truth[idx]
      acc[b, "auroc"]  <- auroc(tr, ext$prob_pos[idx])
      acc[b, "pr_auc"] <- pr_auc(tr, ext$prob_pos[idx])
      acc[b, 3:6]      <- op_metrics(tr, ext$predicted[idx])
    } else {
      rid <- resample_int(); tr <- truth_int[rid]
      if (kind == "leaky") {
        acc[b, "auroc"]  <- mean(vapply(leaky_prob, function(p) auroc(tr, p[rid]), numeric(1)))
        acc[b, "pr_auc"] <- mean(vapply(leaky_prob, function(p) pr_auc(tr, p[rid]), numeric(1)))
        # pooled operating-point metrics across repeats
        prd <- unlist(lapply(reps, function(r) {
          s <- leaky[leaky$repeat_id == r, ]
          setNames(s$predicted, s$sample_id)[rid]
        }))
        acc[b, 3:6] <- op_metrics(rep(tr, length(reps)), prd)
      } else {
        acc[b, "auroc"]  <- auroc(tr, nested_prob[rid])
        acc[b, "pr_auc"] <- pr_auc(tr, nested_prob[rid])
        acc[b, 3:6]      <- op_metrics(tr, nested_pred[rid])
      }
    }
  }
  acc
}

points <- list(leaky = leaky_point, nested = nested_point, external = ext_point)
labels <- c(leaky = "A_leaky_baseline", nested = "B_guarded_nested",
            external = "B_guarded_nested_external")
cohorts <- c(leaky = "GSE25055", nested = "GSE25055", external = "GSE25065")

ci_rows <- list()
for (kind in c("leaky", "nested", "external")) {
  acc <- boot_cohort(kind)
  for (m in METRICS) {
    q <- pct_ci(acc[, m])
    ci_rows[[length(ci_rows) + 1]] <- data.frame(
      cohort = cohorts[[kind]], pipeline = labels[[kind]], metric = m,
      point = round(points[[kind]][[m]], 4),
      ci_lo = round(q[1], 4), ci_hi = round(q[2], 4),
      n_boot = B, method = "stratified percentile bootstrap",
      stringsAsFactors = FALSE)
  }
}
ci_df <- do.call(rbind, ci_rows)
write.csv(ci_df, "tables/uncertainty/bootstrap_ci.csv", row.names = FALSE)

# ---- paired bootstrap for the leakage gap (leaky - guarded) -----------------
d_au <- numeric(B); d_pr <- numeric(B)
for (b in seq_len(B)) {
  rid <- resample_int(); tr <- truth_int[rid]
  la <- mean(vapply(leaky_prob, function(p) auroc(tr, p[rid]), numeric(1)))
  lp <- mean(vapply(leaky_prob, function(p) pr_auc(tr, p[rid]), numeric(1)))
  d_au[b] <- la - auroc(tr, nested_prob[rid])
  d_pr[b] <- lp - pr_auc(tr, nested_prob[rid])
}
boot_p <- function(x) round(2 * min(mean(x > 0), mean(x < 0)), 4)
delta_df <- data.frame(
  comparison = "leaky - guarded (GSE25055)",
  metric = c("delta_auroc", "delta_pr_auc"),
  point = round(c(leaky_point["auroc"] - nested_point["auroc"],
                  leaky_point["pr_auc"] - nested_point["pr_auc"]), 4),
  ci_lo = round(c(pct_ci(d_au)[1], pct_ci(d_pr)[1]), 4),
  ci_hi = round(c(pct_ci(d_au)[2], pct_ci(d_pr)[2]), 4),
  boot_p = c(boot_p(d_au), boot_p(d_pr)),
  n_boot = B, method = "paired stratified bootstrap",
  stringsAsFactors = FALSE)
write.csv(delta_df, "tables/uncertainty/delta_auroc_prauc_ci.csv", row.names = FALSE)

# ---- DeLong test per repeat (leaky_r vs guarded nested, same samples) -------
delong_rows <- lapply(reps, function(r) {
  r1 <- pROC::roc(response = truth_int, predictor = leaky_prob[[r]],
                  levels = c(NEG, POS), direction = "<", quiet = TRUE)
  r2 <- pROC::roc(response = truth_int, predictor = nested_prob,
                  levels = c(NEG, POS), direction = "<", quiet = TRUE)
  tt <- pROC::roc.test(r1, r2, method = "delong", paired = TRUE)
  data.frame(comparison = sprintf("leaky_repeat_%d vs guarded_nested", r),
             auroc_leaky = round(as.numeric(pROC::auc(r1)), 4),
             auroc_guarded = round(as.numeric(pROC::auc(r2)), 4),
             auroc_diff = round(as.numeric(pROC::auc(r1) - pROC::auc(r2)), 4),
             delong_z = round(unname(tt$statistic), 4),
             delong_p = round(tt$p.value, 4), n = length(ids),
             stringsAsFactors = FALSE)
})
delong_df <- do.call(rbind, delong_rows)
write.csv(delong_df, "results/uncertainty/delong_tests.csv", row.names = FALSE)

message("[07] Bootstrap CIs and DeLong tests written.")
message(sprintf("[07] Leakage gap AUROC=%.4f (95%% CI %.4f, %.4f; boot p=%.3f)",
                delta_df$point[1], delta_df$ci_lo[1], delta_df$ci_hi[1], delta_df$boot_p[1]))
message(sprintf("[07] DeLong median p (leaky vs guarded) = %.4f", median(delong_df$delong_p)))
invisible(list(ci = ci_df, delta = delta_df, delong = delong_df))
