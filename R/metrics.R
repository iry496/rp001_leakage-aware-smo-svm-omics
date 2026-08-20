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

# Confusion-matrix metrics under an explicit, predeclared operating point.
# Keeping thresholding outside the model's package-specific class prediction
# makes operating-point comparisons reproducible across cohorts.
compute_operating_metrics <- function(truth, probability, positive,
                                      threshold = 0.5) {
  truth <- as.character(truth)
  negative <- setdiff(unique(truth), positive)
  if (length(negative) != 1L) {
    stop("Operating-point metrics require exactly one negative class.")
  }
  estimate <- ifelse(probability >= threshold, positive, negative)
  tp <- sum(truth == positive & estimate == positive)
  fn <- sum(truth == positive & estimate != positive)
  tn <- sum(truth != positive & estimate != positive)
  fp <- sum(truth != positive & estimate == positive)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  data.frame(
    threshold = threshold,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    mcc = mcc_binary(truth, estimate, positive),
    sensitivity = sensitivity,
    specificity = specificity,
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn,
    stringsAsFactors = FALSE
  )
}

# Probability-quality and calibration diagnostics. The joint calibration
# intercept/slope come from logit(Y) = a + b*logit(p). Calibration-in-the-large
# is estimated separately with logit(p) as an offset. These are descriptive
# external-validation diagnostics, not post-hoc recalibration procedures.
compute_probability_metrics <- function(truth, probability, positive,
                                        clip = 1e-6) {
  truth <- as.character(truth)
  if (length(truth) != length(probability)) {
    stop("truth and probability must have the same length.")
  }
  if (any(!is.finite(probability))) {
    stop("probability contains non-finite values.")
  }
  if (any(probability < 0 | probability > 1)) {
    stop("probability must lie in [0, 1].")
  }

  y <- as.integer(truth == positive)
  p <- pmin(pmax(as.numeric(probability), clip), 1 - clip)
  prevalence <- mean(y)
  brier <- mean((y - p)^2)
  null_brier <- prevalence * (1 - prevalence)
  brier_skill <- if (null_brier > 0) 1 - brier / null_brier else NA_real_
  log_loss <- -mean(y * log(p) + (1 - y) * log(1 - p))
  lp <- stats::qlogis(p)

  joint <- try(suppressWarnings(stats::glm(y ~ lp, family = stats::binomial())),
               silent = TRUE)
  joint_coef <- if (inherits(joint, "try-error")) c(NA_real_, NA_real_) else {
    co <- stats::coef(joint)
    if (length(co) == 2L && all(is.finite(co))) unname(co) else c(NA_real_, NA_real_)
  }
  citl_fit <- try(suppressWarnings(stats::glm(
    y ~ 1 + offset(lp), family = stats::binomial()
  )), silent = TRUE)
  citl <- if (inherits(citl_fit, "try-error")) NA_real_ else {
    co <- stats::coef(citl_fit)
    if (length(co) && is.finite(co[[1]])) unname(co[[1]]) else NA_real_
  }

  data.frame(
    prevalence = prevalence,
    mean_predicted_probability = mean(p),
    expected_observed_ratio = if (prevalence > 0) mean(p) / prevalence else NA_real_,
    brier_score = brier,
    brier_skill_score = brier_skill,
    log_loss = log_loss,
    calibration_intercept = joint_coef[[1]],
    calibration_slope = joint_coef[[2]],
    calibration_in_the_large = citl,
    stringsAsFactors = FALSE
  )
}

# One-row evaluation contract used by the matched ablation and the unified
# external-validation analysis.
evaluate_binary_predictions <- function(truth, probability, positive,
                                        threshold = 0.5) {
  truth <- as.character(truth)
  negative <- setdiff(unique(truth), positive)
  if (length(negative) != 1L) {
    stop("Binary evaluation requires exactly one negative class.")
  }
  estimate <- ifelse(probability >= threshold, positive, negative)
  discrim <- compute_binary_metrics(
    truth = factor(truth, levels = c(negative, positive)),
    estimate = factor(estimate, levels = c(negative, positive)),
    probability = probability,
    positive = positive
  )
  operating <- compute_operating_metrics(truth, probability, positive, threshold)
  probability_metrics <- compute_probability_metrics(truth, probability, positive)
  data.frame(
    n = discrim$n,
    auroc = discrim$auroc,
    pr_auc = discrim$pr_auc,
    operating[, c("threshold", "balanced_accuracy", "mcc", "sensitivity",
                  "specificity", "tp", "fp", "tn", "fn")],
    probability_metrics,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
