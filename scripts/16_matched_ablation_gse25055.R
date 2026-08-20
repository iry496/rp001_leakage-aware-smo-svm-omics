#!/usr/bin/env Rscript
# =============================================================================
# Matched one-factor ablation: supervised feature-selection placement only.
#
# Both arms use identical samples, repeated outer folds, inner folds, top-K,
# cost grid, tuning criterion, scaling, class weights, probability procedure,
# model family, threshold, and performance estimators. The sole design change:
#
#   leaky   : t-test top-K selected once on all GSE25055 samples before CV
#   guarded : t-test top-K refit only inside the relevant training partition
#
# Run:
#   Rscript scripts/16_matched_ablation_gse25055.R smoke 3
#   Rscript scripts/16_matched_ablation_gse25055.R full 30
# =============================================================================

suppressWarnings(suppressMessages({
  source("R/00_config.R")
  source("R/feature_selection.R")
  source("R/preprocessing.R")
  source("R/model_smo_svm.R")
  source("R/metrics.R")
  source("R/leakage_checks.R")
}))

args <- commandArgs(trailingOnly = TRUE)
RUN_TAG <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else "smoke"
N_REPEATS <- if (length(args) >= 2L) as.integer(args[[2]]) else 3L
if (!is.finite(N_REPEATS) || N_REPEATS < 1L) stop("N_REPEATS must be positive.")

ACCESSION <- "GSE25055"
LABEL_FIELD <- "pathologic_response_pcr_rd"
POSITIVE <- "pCR"
NEGATIVE <- "RD"
TOP_K <- 100L
OUTER_FOLDS <- 5L
INNER_FOLDS <- 5L
COST_GRID <- c(0.25, 1, 4)
KERNEL <- "linear"
THRESHOLD <- 0.5
BOOTSTRAP_B <- 2000L
N_WORKERS <- max(1L, min(4L, parallel::detectCores() - 2L))

RESULTS_DIR <- file.path("results", "matched_ablation")
TABLES_DIR <- file.path("tables", "matched_ablation")
FIGURES_DIR <- file.path("figures", "matched_ablation")
CACHE_X <- file.path("processed_data", "gse25055_perm_cache.rds")
CACHE_Y <- file.path("processed_data", "gse25055_perm_labels.rds")
for (d in c(RESULTS_DIR, TABLES_DIR, FIGURES_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

load_gse25055 <- function() {
  if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery required.")
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase required.")
  g <- GEOquery::getGEO(ACCESSION, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  expr <- Biobase::exprs(g)
  pheno <- Biobase::pData(g)
  label_col <- NULL
  for (cn in colnames(pheno)) {
    values <- as.character(pheno[[cn]])
    if (any(grepl(LABEL_FIELD, values, ignore.case = TRUE)) ||
        grepl(LABEL_FIELD, cn, ignore.case = TRUE)) {
      label_col <- cn
      break
    }
  }
  if (is.null(label_col)) stop("GSE25055 endpoint field not found.")
  raw <- trimws(sub(
    paste0("^.*", LABEL_FIELD, "\\s*[:=]?\\s*"), "",
    as.character(pheno[[label_col]]), ignore.case = TRUE
  ))
  keep <- !(is.na(raw) | raw %in% c("NA", "na", "N/A", "", "NaN"))
  y <- factor(raw[keep], levels = c(NEGATIVE, POSITIVE))
  if (any(is.na(y))) stop("Unexpected GSE25055 endpoint coding.")
  x <- t(expr[, keep, drop = FALSE])
  rownames(x) <- colnames(expr)[keep]
  ord <- order(rownames(x))
  list(x = x[ord, , drop = FALSE], y = y[ord])
}

load_or_cache <- function() {
  if (file.exists(CACHE_X) && file.exists(CACHE_Y)) {
    x <- readRDS(CACHE_X)
    y <- readRDS(CACHE_Y)
    ord <- order(rownames(x))
    return(list(x = x[ord, , drop = FALSE], y = y[ord]))
  }
  d <- load_gse25055()
  x <- filter_near_zero_variance(d$x, cutoff = 1e-8)$x
  dir.create("processed_data", recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, CACHE_X)
  saveRDS(d$y, CACHE_Y)
  list(x = x, y = d$y)
}

# Vectorized Welch selector; validated at runtime against the canonical t.test
# implementation before the ablation is allowed to run.
select_ttest_fast <- function(x, y, top_k = TOP_K) {
  x <- as.matrix(x)
  y <- as.factor(y)
  lev <- levels(y)
  g1 <- which(y == lev[[1]])
  g2 <- which(y == lev[[2]])
  n1 <- length(g1)
  n2 <- length(g2)
  m1 <- colMeans(x[g1, , drop = FALSE])
  m2 <- colMeans(x[g2, , drop = FALSE])
  v1 <- (colSums(x[g1, , drop = FALSE]^2) - n1 * m1^2) / (n1 - 1)
  v2 <- (colSums(x[g2, , drop = FALSE]^2) - n2 * m2^2) / (n2 - 1)
  v1[v1 < 0] <- 0
  v2[v2 < 0] <- 0
  score <- abs((m1 - m2) / sqrt(v1 / n1 + v2 / n2))
  score[!is.finite(score)] <- 0
  names(score) <- colnames(x)
  names(sort(score, decreasing = TRUE))[seq_len(min(top_k, length(score)))]
}

validate_selector <- function(x, y) {
  cases <- c(list(observed = y), setNames(lapply(seq_len(5L), function(i) {
    set.seed(SEED + 900000L + i)
    y[sample.int(length(y))]
  }), paste0("permuted_", seq_len(5L))))
  out <- do.call(rbind, lapply(names(cases), function(label) {
    yy <- cases[[label]]
    canonical <- select_features(x, yy, method = "t_test", top_k = TOP_K)$features
    fast <- select_ttest_fast(x, yy, TOP_K)
    data.frame(
      case = label,
      overlap = length(intersect(canonical, fast)),
      exact_order = identical(canonical, fast),
      stringsAsFactors = FALSE
    )
  }))
  if (!all(out$exact_order)) stop("Fast selector failed canonical equivalence gate.")
  out
}

balanced_accuracy_at_threshold <- function(truth, probability) {
  compute_operating_metrics(truth, probability, POSITIVE, THRESHOLD)$balanced_accuracy
}

make_outer_folds <- function(y, repeat_id) {
  set.seed(SEED + repeat_id - 1L)
  caret::createFolds(y, k = OUTER_FOLDS, returnTrain = FALSE)
}

make_inner_folds <- function(y_train, repeat_id, outer_fold) {
  set.seed(SEED + 100000L + repeat_id * 100L + outer_fold)
  caret::createFolds(y_train, k = INNER_FOLDS, returnTrain = FALSE)
}

model_seed <- function(repeat_id, outer_fold, cost_index, inner_fold = 0L) {
  SEED + 200000L + repeat_id * 10000L + outer_fold * 1000L +
    cost_index * 100L + inner_fold
}

final_model_seed <- function(repeat_id, outer_fold) {
  # Identical across arms even when the common tuning policy selects different
  # costs downstream of the feature-placement intervention.
  SEED + 700000L + repeat_id * 100L + outer_fold
}

tune_cost <- function(x_train, y_train, inner_folds, arm, global_features,
                      repeat_id, outer_fold) {
  scores <- numeric(length(COST_GRID))
  for (ci in seq_along(COST_GRID)) {
    fold_scores <- numeric(length(inner_folds))
    for (inner_fold in seq_along(inner_folds)) {
      validation <- inner_folds[[inner_fold]]
      training <- setdiff(seq_along(y_train), validation)
      features <- if (identical(arm, "leaky_global_selection")) {
        global_features
      } else {
        select_ttest_fast(
          x_train[training, , drop = FALSE], y_train[training], TOP_K
        )
      }
      scaler <- fit_standardizer(x_train[training, features, drop = FALSE])
      x_tr <- apply_standardizer(x_train[training, features, drop = FALSE], scaler)
      x_va <- apply_standardizer(x_train[validation, features, drop = FALSE], scaler)
      set.seed(model_seed(repeat_id, outer_fold, ci, inner_fold))
      model <- train_smo_svm(
        x_tr, y_train[training], kernel = KERNEL, cost = COST_GRID[[ci]],
        class_weights = TRUE, probability = TRUE
      )
      probability <- predict_smo_svm(
        model, x_va, positive_class = POSITIVE
      )$probability_positive
      if (is.null(probability)) stop("SVM probability estimation failed.")
      fold_scores[[inner_fold]] <- balanced_accuracy_at_threshold(
        y_train[validation], probability
      )
    }
    scores[[ci]] <- mean(fold_scores, na.rm = TRUE)
  }
  list(
    best_cost = COST_GRID[[which.max(scores)]],
    scores = data.frame(cost = COST_GRID, mean_inner_balanced_accuracy = scores)
  )
}

fit_outer_arm <- function(x, y, train_idx, test_idx, inner_folds, arm,
                          global_features, repeat_id, outer_fold) {
  x_train <- x[train_idx, , drop = FALSE]
  y_train <- y[train_idx]
  tuning <- tune_cost(
    x_train, y_train, inner_folds, arm, global_features, repeat_id, outer_fold
  )
  features <- if (identical(arm, "leaky_global_selection")) {
    global_features
  } else {
    select_ttest_fast(x_train, y_train, TOP_K)
  }
  scaler <- fit_standardizer(x_train[, features, drop = FALSE])
  x_tr <- apply_standardizer(x_train[, features, drop = FALSE], scaler)
  x_te <- apply_standardizer(x[test_idx, features, drop = FALSE], scaler)
  set.seed(final_model_seed(repeat_id, outer_fold))
  model <- train_smo_svm(
    x_tr, y_train, kernel = KERNEL, cost = tuning$best_cost,
    class_weights = TRUE, probability = TRUE
  )
  probability <- predict_smo_svm(
    model, x_te, positive_class = POSITIVE
  )$probability_positive
  prediction <- ifelse(probability >= THRESHOLD, POSITIVE, NEGATIVE)
  list(
    predictions = data.frame(
      repeat_id = repeat_id,
      outer_seed = SEED + repeat_id - 1L,
      outer_fold = outer_fold,
      arm = arm,
      sample_id = rownames(x)[test_idx],
      truth = as.character(y[test_idx]),
      probability = probability,
      predicted = prediction,
      stringsAsFactors = FALSE
    ),
    features = data.frame(
      repeat_id = repeat_id,
      outer_fold = outer_fold,
      arm = arm,
      rank = seq_along(features),
      feature = features,
      stringsAsFactors = FALSE
    ),
    tuning = data.frame(
      repeat_id = repeat_id,
      outer_fold = outer_fold,
      arm = arm,
      best_cost = tuning$best_cost,
      tuning$scores,
      stringsAsFactors = FALSE
    )
  )
}

run_repeat <- function(repeat_id, x, y, global_features) {
  start <- Sys.time()
  folds <- make_outer_folds(y, repeat_id)
  pieces <- list()
  assignments <- list()
  for (outer_fold in seq_along(folds)) {
    test_idx <- folds[[outer_fold]]
    train_idx <- setdiff(seq_along(y), test_idx)
    assert_no_overlap(rownames(x)[train_idx], rownames(x)[test_idx])
    inner_folds <- make_inner_folds(y[train_idx], repeat_id, outer_fold)
    assignments[[outer_fold]] <- data.frame(
      repeat_id = repeat_id,
      outer_seed = SEED + repeat_id - 1L,
      outer_fold = outer_fold,
      sample_id = rownames(x)[test_idx],
      stringsAsFactors = FALSE
    )
    for (arm in c("leaky_global_selection", "guarded_within_partition")) {
      pieces[[length(pieces) + 1L]] <- fit_outer_arm(
        x, y, train_idx, test_idx, inner_folds, arm, global_features,
        repeat_id, outer_fold
      )
    }
  }
  list(
    predictions = do.call(rbind, lapply(pieces, `[[`, "predictions")),
    features = do.call(rbind, lapply(pieces, `[[`, "features")),
    tuning = do.call(rbind, lapply(pieces, `[[`, "tuning")),
    assignments = do.call(rbind, assignments),
    runtime_seconds = as.numeric(difftime(Sys.time(), start, units = "secs"))
  )
}

metric_names <- c(
  "auroc", "pr_auc", "balanced_accuracy", "mcc", "sensitivity",
  "specificity", "brier_score"
)

evaluate_group <- function(d) {
  evaluate_binary_predictions(d$truth, d$probability, POSITIVE, THRESHOLD)
}

per_repeat_metrics <- function(predictions) {
  groups <- split(predictions, interaction(predictions$repeat_id, predictions$arm,
                                           drop = TRUE))
  long <- do.call(rbind, lapply(groups, function(d) {
    m <- evaluate_group(d)
    data.frame(
      repeat_id = d$repeat_id[[1]], arm = d$arm[[1]],
      m[, metric_names, drop = FALSE], stringsAsFactors = FALSE
    )
  }))
  rownames(long) <- NULL
  leaky <- long[long$arm == "leaky_global_selection", ]
  guarded <- long[long$arm == "guarded_within_partition", ]
  leaky <- leaky[order(leaky$repeat_id), ]
  guarded <- guarded[order(guarded$repeat_id), ]
  if (!identical(leaky$repeat_id, guarded$repeat_id)) stop("Unmatched repeat metrics.")
  wide <- data.frame(repeat_id = leaky$repeat_id)
  for (metric in metric_names) {
    wide[[paste0("leaky_", metric)]] <- leaky[[metric]]
    wide[[paste0("guarded_", metric)]] <- guarded[[metric]]
    wide[[paste0("delta_", metric)]] <- leaky[[metric]] - guarded[[metric]]
  }
  list(long = long, wide = wide)
}

summarize_paired_differences <- function(wide) {
  do.call(rbind, lapply(metric_names, function(metric) {
    delta <- wide[[paste0("delta_", metric)]]
    data.frame(
      metric = metric,
      n_repeats = length(delta),
      leaky_mean = mean(wide[[paste0("leaky_", metric)]], na.rm = TRUE),
      guarded_mean = mean(wide[[paste0("guarded_", metric)]], na.rm = TRUE),
      mean_delta_leaky_minus_guarded = mean(delta, na.rm = TRUE),
      median_delta = median(delta, na.rm = TRUE),
      sd_delta = stats::sd(delta, na.rm = TRUE),
      empirical_q025 = unname(stats::quantile(delta, 0.025, na.rm = TRUE)),
      empirical_q975 = unname(stats::quantile(delta, 0.975, na.rm = TRUE)),
      n_delta_positive = sum(delta > 0, na.rm = TRUE),
      interval_label = "empirical split-distribution interval; not an independence-based CI",
      stringsAsFactors = FALSE
    )
  }))
}

aggregate_cross_fitted_predictions <- function(predictions) {
  key <- interaction(predictions$arm, predictions$sample_id, drop = TRUE)
  out <- do.call(rbind, lapply(split(predictions, key), function(d) {
    data.frame(
      arm = d$arm[[1]], sample_id = d$sample_id[[1]], truth = d$truth[[1]],
      mean_probability = mean(d$probability),
      n_repeats = length(unique(d$repeat_id)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

paired_patient_bootstrap <- function(aggregated, B = BOOTSTRAP_B) {
  leaky <- aggregated[aggregated$arm == "leaky_global_selection", ]
  guarded <- aggregated[aggregated$arm == "guarded_within_partition", ]
  leaky <- leaky[order(leaky$sample_id), ]
  guarded <- guarded[order(guarded$sample_id), ]
  if (!identical(leaky$sample_id, guarded$sample_id) ||
      !identical(leaky$truth, guarded$truth)) stop("Patient-level arm alignment failed.")
  truth <- leaky$truth
  pos <- which(truth == POSITIVE)
  neg <- which(truth == NEGATIVE)
  point_l <- evaluate_binary_predictions(truth, leaky$mean_probability, POSITIVE, THRESHOLD)
  point_g <- evaluate_binary_predictions(truth, guarded$mean_probability, POSITIVE, THRESHOLD)
  set.seed(SEED + 800000L)
  boot <- array(NA_real_, dim = c(B, length(metric_names), 3L),
                dimnames = list(NULL, metric_names, c("leaky", "guarded", "delta")))
  for (b in seq_len(B)) {
    idx <- c(sample(pos, replace = TRUE), sample(neg, replace = TRUE))
    ml <- evaluate_binary_predictions(
      truth[idx], leaky$mean_probability[idx], POSITIVE, THRESHOLD
    )
    mg <- evaluate_binary_predictions(
      truth[idx], guarded$mean_probability[idx], POSITIVE, THRESHOLD
    )
    boot[b, , "leaky"] <- as.numeric(ml[1, metric_names])
    boot[b, , "guarded"] <- as.numeric(mg[1, metric_names])
    boot[b, , "delta"] <- boot[b, , "leaky"] - boot[b, , "guarded"]
  }
  do.call(rbind, lapply(metric_names, function(metric) {
    ci <- stats::quantile(boot[, metric, "delta"], c(0.025, 0.975), na.rm = TRUE)
    data.frame(
      metric = metric,
      leaky_point = point_l[[metric]],
      guarded_point = point_g[[metric]],
      delta_leaky_minus_guarded = point_l[[metric]] - point_g[[metric]],
      delta_ci_lo = unname(ci[[1]]),
      delta_ci_hi = unname(ci[[2]]),
      n_boot = B,
      method = paste(
        "paired stratified patient bootstrap of repeat-averaged cross-fitted",
        "predictions; conditional on realized CV fits; no workflow refitting"
      ),
      stringsAsFactors = FALSE
    )
  }))
}

write_figure <- function(summary, paired) {
  auroc <- summary[summary$metric == "auroc", ]
  pr <- summary[summary$metric == "pr_auc", ]
  values <- c(auroc$mean_delta_leaky_minus_guarded, pr$mean_delta_leaky_minus_guarded)
  lows <- c(auroc$empirical_q025, pr$empirical_q025)
  highs <- c(auroc$empirical_q975, pr$empirical_q975)
  for (device in list(
    function() grDevices::png(file.path(FIGURES_DIR, paste0("matched_ablation_", RUN_TAG, ".png")),
                              width = 1400, height = 650, res = 120),
    function() grDevices::pdf(file.path(FIGURES_DIR, paste0("matched_ablation_", RUN_TAG, ".pdf")),
                              width = 11.5, height = 5.3)
  )) {
    device()
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    ylim <- range(c(paired$leaky_auroc, paired$guarded_auroc))
    graphics::plot(c(1, 2), c(NA, NA), xlim = c(0.75, 2.25), ylim = ylim,
                   xaxt = "n", xlab = "", ylab = "AUROC",
                   main = "Thirty paired repeated-CV estimates")
    for (i in seq_len(nrow(paired))) {
      graphics::segments(1, paired$leaky_auroc[[i]], 2, paired$guarded_auroc[[i]],
                         col = grDevices::adjustcolor("#64748B", alpha.f = 0.32))
    }
    graphics::points(rep(1, nrow(paired)), paired$leaky_auroc,
                     pch = 19, col = grDevices::adjustcolor("#B91C1C", alpha.f = 0.55))
    graphics::points(rep(2, nrow(paired)), paired$guarded_auroc,
                     pch = 19, col = grDevices::adjustcolor("#047857", alpha.f = 0.55))
    graphics::points(c(1, 2), c(mean(paired$leaky_auroc), mean(paired$guarded_auroc)),
                     pch = 21, cex = 1.7, lwd = 2, bg = c("#FCA5A5", "#6EE7B7"))
    graphics::axis(1, at = c(1, 2), labels = c("Global selection\n(leaky)",
                                               "Within-partition\n(guarded)"))

    set.seed(SEED)
    graphics::plot(seq_along(values), values, pch = 19, cex = 1.2, xaxt = "n",
                   ylim = range(c(lows, highs, 0)), xlab = "",
                   ylab = "Paired difference (leaky - guarded)",
                   main = "Effect of feature-selection placement")
    graphics::points(jitter(rep(1, nrow(paired)), amount = 0.055), paired$delta_auroc,
                     pch = 16, cex = 0.65,
                     col = grDevices::adjustcolor("#1D4ED8", alpha.f = 0.45))
    graphics::points(jitter(rep(2, nrow(paired)), amount = 0.055), paired$delta_pr_auc,
                     pch = 16, cex = 0.65,
                     col = grDevices::adjustcolor("#7C3AED", alpha.f = 0.45))
    graphics::arrows(seq_along(values), lows, seq_along(values), highs,
                     angle = 90, code = 3, length = 0.06, lwd = 2)
    graphics::points(seq_along(values), values, pch = 21, cex = 1.5,
                     lwd = 2, bg = "white")
    graphics::axis(1, at = seq_along(values), labels = c("AUROC", "PR-AUC"))
    graphics::abline(h = 0, lty = 3)
    graphics::mtext("Bars: empirical 2.5th-97.5th percentile split interval",
                    side = 1, line = 3, cex = 0.78)
    grDevices::dev.off()
  }
}

main <- function() {
  start <- Sys.time()
  if (!requireNamespace("caret", quietly = TRUE)) stop("caret required.")
  d <- load_or_cache()
  x <- d$x
  y <- d$y
  selector_check <- validate_selector(x, y)
  global_features <- select_ttest_fast(x, y, TOP_K)
  repeat_ids <- seq_len(N_REPEATS)
  message(sprintf(
    "[matched] n=%d p=%d repeats=%d workers=%d", nrow(x), ncol(x),
    N_REPEATS, N_WORKERS
  ))
  worker <- function(i) run_repeat(i, x, y, global_features)
  if (.Platform$OS.type != "windows" && N_WORKERS > 1L) {
    pieces <- parallel::mclapply(
      repeat_ids, worker, mc.cores = N_WORKERS, mc.preschedule = FALSE
    )
  } else {
    pieces <- lapply(repeat_ids, worker)
  }
  predictions <- do.call(rbind, lapply(pieces, `[[`, "predictions"))
  features <- do.call(rbind, lapply(pieces, `[[`, "features"))
  tuning <- do.call(rbind, lapply(pieces, `[[`, "tuning"))
  assignments <- do.call(rbind, lapply(pieces, `[[`, "assignments"))

  # Hard gate: the held-out sample/fold keys must match exactly between arms.
  lk <- predictions[predictions$arm == "leaky_global_selection",
                    c("repeat_id", "outer_fold", "sample_id", "truth")]
  gd <- predictions[predictions$arm == "guarded_within_partition",
                    c("repeat_id", "outer_fold", "sample_id", "truth")]
  lk <- lk[do.call(order, lk), ]
  gd <- gd[do.call(order, gd), ]
  rownames(lk) <- NULL
  rownames(gd) <- NULL
  if (!identical(lk, gd)) stop("Matched outer-fold identity gate failed.")

  metrics <- per_repeat_metrics(predictions)
  paired_summary <- summarize_paired_differences(metrics$wide)
  aggregated <- aggregate_cross_fitted_predictions(predictions)
  bootstrap <- paired_patient_bootstrap(aggregated)
  write_figure(paired_summary, metrics$wide)

  contract <- data.frame(
    component = c(
      "samples", "outer_folds", "inner_folds", "feature_budget", "selector",
      "cost_grid", "tuning_metric", "scaling", "class_weights", "model",
      "probability_procedure", "final_fit_rng_seed", "threshold", "performance_estimators",
      "supervised_feature_selection_placement"
    ),
    leaky = c(
      "GSE25055 n=306", "shared", "shared", TOP_K, "absolute Welch t statistic",
      paste(COST_GRID, collapse = "|"), "balanced accuracy at 0.5",
      "training-partition center/scale", "inverse frequency", "linear e1071 SVM",
      "libsvm probability=TRUE", "shared by repeat and outer fold", THRESHOLD, "shared",
      "once on all samples before CV"
    ),
    guarded = c(
      "GSE25055 n=306", "shared", "shared", TOP_K, "absolute Welch t statistic",
      paste(COST_GRID, collapse = "|"), "balanced accuracy at 0.5",
      "training-partition center/scale", "inverse frequency", "linear e1071 SVM",
      "libsvm probability=TRUE", "shared by repeat and outer fold", THRESHOLD, "shared",
      "refit inside each training partition"
    ),
    matched = c(rep(TRUE, 14L), FALSE),
    stringsAsFactors = FALSE
  )

  write.csv(predictions, file.path(RESULTS_DIR, paste0("predictions_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(features, file.path(RESULTS_DIR, paste0("selected_features_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(tuning, file.path(RESULTS_DIR, paste0("tuning_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(assignments, file.path(RESULTS_DIR, paste0("fold_assignments_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(selector_check, file.path(RESULTS_DIR, "selector_equivalence_check.csv"), row.names = FALSE)
  write.csv(contract, file.path(TABLES_DIR, "matched_design_contract.csv"), row.names = FALSE)
  write.csv(metrics$long, file.path(TABLES_DIR, paste0("per_repeat_metrics_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(metrics$wide, file.path(TABLES_DIR, paste0("paired_differences_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(paired_summary, file.path(TABLES_DIR, paste0("paired_summary_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(aggregated, file.path(RESULTS_DIR, paste0("repeat_averaged_predictions_", RUN_TAG, ".csv")), row.names = FALSE)
  write.csv(bootstrap, file.path(TABLES_DIR, paste0("paired_patient_bootstrap_", RUN_TAG, ".csv")), row.names = FALSE)

  a <- paired_summary[paired_summary$metric == "auroc", ]
  p <- paired_summary[paired_summary$metric == "pr_auc", ]
  notes <- c(
    paste0("# Matched one-factor ablation - ", RUN_TAG), "",
    sprintf("- GSE25055: n=%d; repeated nested CV: %d x %d outer folds, %d inner folds.",
            nrow(x), N_REPEATS, OUTER_FOLDS, INNER_FOLDS),
    "- The only unmatched component is supervised feature-selection placement.",
    "- All outer test-fold keys passed exact arm-identity checks.",
    "- Final-model RNG seeds were identical across arms within every repeat and outer fold.",
    "- Fast Welch selector passed exact-order equivalence against the canonical selector.",
    "", "## Paired split-level results",
    sprintf("- AUROC: leaky mean %.4f; guarded mean %.4f; mean delta %.4f (empirical split interval %.4f to %.4f; positive in %d/%d repeats).",
            a$leaky_mean, a$guarded_mean, a$mean_delta_leaky_minus_guarded,
            a$empirical_q025, a$empirical_q975, a$n_delta_positive, a$n_repeats),
    sprintf("- PR-AUC: leaky mean %.4f; guarded mean %.4f; mean delta %.4f (empirical split interval %.4f to %.4f; positive in %d/%d repeats).",
            p$leaky_mean, p$guarded_mean, p$mean_delta_leaky_minus_guarded,
            p$empirical_q025, p$empirical_q975, p$n_delta_positive, p$n_repeats),
    "- Split repeats reuse the same patients and are correlated; the empirical interval is a robustness distribution, not an independence-based confidence interval.",
    "- Patient-bootstrap intervals use repeat-averaged cross-fitted scores and are conditional on the realized CV fits; the workflow is not refit within bootstrap samples.",
    "", "## Interpretation",
    "This matched contrast isolates supervised feature-selection placement under the stated pipeline contract. It replaces the older unmatched complete-workflow contrast as the leakage-effect estimate; the older comparison remains useful only as a complete-workflow audit.",
    "", sprintf("_Runtime %.1f minutes._", as.numeric(difftime(Sys.time(), start, units = "mins")))
  )
  writeLines(notes, file.path(RESULTS_DIR, paste0("matched_ablation_", RUN_TAG, "_notes.md")))
  message(sprintf("[matched] complete in %.1f minutes.", as.numeric(difftime(Sys.time(), start, units = "mins"))))
  invisible(list(summary = paired_summary, bootstrap = bootstrap))
}

if (sys.nframe() == 0L) main()
