# Performance metrics for imbalanced binary classification.

mcc_binary <- function(truth, estimate, positive) {
  truth <- factor(truth)
  estimate <- factor(estimate, levels = levels(truth))
  tp <- as.numeric(sum(truth == positive & estimate == positive))
  tn <- as.numeric(sum(truth != positive & estimate != positive))
  fp <- as.numeric(sum(truth != positive & estimate == positive))
  fn <- as.numeric(sum(truth == positive & estimate != positive))
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  if (denom == 0) return(NA_real_)
  (tp * tn - fp * fn) / denom
}

balanced_accuracy_binary <- function(truth, estimate, positive) {
  truth <- factor(truth)
  estimate <- factor(estimate, levels = levels(truth))
  sens <- sum(truth == positive & estimate == positive) / sum(truth == positive)
  spec <- sum(truth != positive & estimate != positive) / sum(truth != positive)
  mean(c(sens, spec), na.rm = TRUE)
}

compute_binary_metrics <- function(truth, estimate, probability = NULL, positive) {
  out <- data.frame(
    n = length(truth),
    mcc = mcc_binary(truth, estimate, positive),
    balanced_accuracy = balanced_accuracy_binary(truth, estimate, positive),
    stringsAsFactors = FALSE
  )
  if (!is.null(probability)) {
    if (requireNamespace("pROC", quietly = TRUE)) {
      # Fix ROC orientation a priori. pROC's automatic direction selection uses
      # the observed data and can bias AUC upward in resampling/permutation
      # analyses. Probability is the model's score for `positive`, so the
      # negative class is the control level and increasing scores indicate the
      # positive class.
      negative <- setdiff(levels(factor(truth)), positive)
      if (length(negative) != 1L) stop("Binary AUROC requires exactly one negative class.")
      out$auroc <- as.numeric(pROC::auc(pROC::roc(
        response = truth,
        predictor = probability,
        levels = c(negative, positive),
        direction = "<",
        quiet = TRUE
      )))
    } else {
      out$auroc <- NA_real_
    }
    if (requireNamespace("PRROC", quietly = TRUE)) {
      truth_binary <- as.integer(truth == positive)
      out$pr_auc <- tryCatch({
        PRROC::pr.curve(scores.class0 = probability[truth_binary == 1],
                        scores.class1 = probability[truth_binary == 0], curve = FALSE)$auc.integral
      }, error = function(e) NA_real_)
    } else {
      out$pr_auc <- NA_real_
    }
  }
  out
}
