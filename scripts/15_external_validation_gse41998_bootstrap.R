#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pROC)
  library(PRROC)
})

set.seed(20260620)
B <- 2000L
POS <- "pCR"
NEG <- "RD"

d <- read.csv(
  "results/external_validation_gse41998/gse41998_predictions.csv",
  stringsAsFactors = FALSE
)

auroc <- function(truth, prob) {
  as.numeric(pROC::auc(pROC::roc(
    response = truth,
    predictor = prob,
    levels = c(NEG, POS),
    direction = "<",
    quiet = TRUE
  )))
}

pr_auc <- function(truth, prob) {
  positive <- truth == POS
  PRROC::pr.curve(
    scores.class0 = prob[positive],
    scores.class1 = prob[!positive],
    curve = FALSE
  )$auc.integral
}

operating_metrics <- function(truth, pred) {
  tp <- as.numeric(sum(truth == POS & pred == POS))
  tn <- as.numeric(sum(truth == NEG & pred == NEG))
  fp <- as.numeric(sum(truth == NEG & pred == POS))
  fn <- as.numeric(sum(truth == POS & pred == NEG))
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  c(
    balanced_accuracy = mean(c(sensitivity, specificity)),
    mcc = if (den == 0) NA_real_ else (tp * tn - fp * fn) / den,
    sensitivity = sensitivity,
    specificity = specificity
  )
}

metric_vector <- function(truth, prob, pred) {
  c(auroc = auroc(truth, prob), pr_auc = pr_auc(truth, prob),
    operating_metrics(truth, pred))
}

pos <- which(d$truth == POS)
neg <- which(d$truth == NEG)
variants <- list(
  discovery_derived_primary = c("prob_pos_primary", "pred_primary"),
  within_cohort_zscore_sensitivity = c("prob_pos_sensitivity", "pred_sensitivity")
)

rows <- list()
for (variant in names(variants)) {
  probability <- variants[[variant]][1]
  prediction <- variants[[variant]][2]
  point <- metric_vector(d$truth, d[[probability]], d[[prediction]])
  boot <- matrix(NA_real_, nrow = B, ncol = length(point),
                 dimnames = list(NULL, names(point)))
  for (b in seq_len(B)) {
    idx <- c(sample(pos, replace = TRUE), sample(neg, replace = TRUE))
    boot[b, ] <- metric_vector(d$truth[idx], d[[probability]][idx], d[[prediction]][idx])
  }
  for (metric in names(point)) {
    ci <- quantile(boot[, metric], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = "GSE41998",
      scaling = variant,
      metric = metric,
      point = round(point[[metric]], 4),
      ci_lo = round(ci[1], 4),
      ci_hi = round(ci[2], 4),
      n_boot = B,
      method = "stratified percentile bootstrap",
      stringsAsFactors = FALSE
    )
  }
}

out <- do.call(rbind, rows)
write.csv(
  out,
  "tables/external_validation_gse41998/gse41998_bootstrap_ci.csv",
  row.names = FALSE
)
message("[15] GSE41998 bootstrap confidence intervals written.")
