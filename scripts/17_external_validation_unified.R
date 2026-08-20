#!/usr/bin/env Rscript
# =============================================================================
# Unified frozen-model external validation.
#
# Fits one GSE25055 discovery model once, freezes it, then applies that exact
# model/probability fit to GSE25065 and GSE41998. Both primary projections use
# discovery-derived scaling and the same explicit P(pCR) >= 0.5 rule.
# GSE41998 also retains the predeclared label-blind within-cohort z-score
# sensitivity analysis. External labels are used only after prediction.
#
# Adds probability-quality/calibration diagnostics and stratified percentile
# bootstrap intervals conditional on the single frozen model and stored scores.
# The bootstrap does not repeat discovery fitting.
# =============================================================================

suppressWarnings(suppressMessages({
  source("R/00_config.R")
  source("R/feature_selection.R")
  source("R/preprocessing.R")
  source("R/model_smo_svm.R")
  source("R/metrics.R")
  source("R/leakage_checks.R")
}))

DISCOVERY <- "GSE25055"
POSITIVE <- "pCR"
NEGATIVE <- "RD"
LABEL_FIELD <- "pathologic_response_pcr_rd"
TOP_K <- 100L
CV_FOLDS <- 5L
COST_GRID <- c(0.25, 1, 4)
KERNEL <- "linear"
THRESHOLD <- 0.5
BOOTSTRAP_B <- 2000L
CACHE_X <- file.path("processed_data", "gse25055_perm_cache.rds")
CACHE_Y <- file.path("processed_data", "gse25055_perm_labels.rds")

RESULTS_DIR <- file.path("results", "external_validation_unified")
TABLES_DIR <- file.path("tables", "external_validation_unified")
FIGURES_DIR <- file.path("figures", "external_validation_unified")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

load_binary_gse <- function(accession, label_field = LABEL_FIELD) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery required.")
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase required.")
  g <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  expr <- Biobase::exprs(g)
  pheno <- Biobase::pData(g)
  label_col <- NULL
  for (cn in colnames(pheno)) {
    values <- as.character(pheno[[cn]])
    if (any(grepl(label_field, values, ignore.case = TRUE)) ||
        grepl(label_field, cn, ignore.case = TRUE)) {
      label_col <- cn
      break
    }
  }
  if (is.null(label_col)) stop("Endpoint field not found in ", accession, ".")
  raw <- trimws(sub(
    paste0("^.*", label_field, "\\s*[:=]?\\s*"), "",
    as.character(pheno[[label_col]]), ignore.case = TRUE
  ))
  keep <- !(is.na(raw) | raw %in% c("NA", "na", "N/A", "", "NaN"))
  y <- factor(raw[keep], levels = c(NEGATIVE, POSITIVE))
  if (any(is.na(y))) stop("Unexpected endpoint coding in ", accession, ".")
  x <- t(expr[, keep, drop = FALSE])
  rownames(x) <- colnames(expr)[keep]
  ord <- order(rownames(x))
  list(
    x = x[ord, , drop = FALSE], y = y[ord], label_col = label_col,
    total_n = length(raw), excluded_n = sum(!keep)
  )
}

load_discovery <- function() {
  if (file.exists(CACHE_X) && file.exists(CACHE_Y)) {
    x <- readRDS(CACHE_X)
    y <- readRDS(CACHE_Y)
    ord <- order(rownames(x))
    return(list(x = x[ord, , drop = FALSE], y = y[ord]))
  }
  d <- load_binary_gse(DISCOVERY)
  x <- filter_near_zero_variance(d$x, cutoff = 1e-8)$x
  dir.create("processed_data", recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, CACHE_X)
  saveRDS(d$y, CACHE_Y)
  list(x = x, y = d$y)
}

load_gse41998 <- function() {
  accession <- "GSE41998"
  if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery required.")
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase required.")
  g <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  expr <- Biobase::exprs(g)
  pheno <- Biobase::pData(g)
  is_pcr_col <- function(cn) {
    v <- tolower(trimws(as.character(pheno[[cn]])))
    mean(grepl("^pcr\\s*:", v)) > 0.5
  }
  candidates <- colnames(pheno)[vapply(colnames(pheno), is_pcr_col, logical(1))]
  if (!length(candidates)) stop("GSE41998 exact pcr: endpoint field not found.")
  label_col <- candidates[[1]]
  bare <- tolower(trimws(sub(
    "^\\s*pcr\\s*:\\s*", "", as.character(pheno[[label_col]]),
    ignore.case = TRUE
  )))
  mapped <- ifelse(bare %in% c("yes", "y"), POSITIVE,
                   ifelse(bare %in% c("no", "n"), NEGATIVE, NA_character_))
  keep <- !is.na(mapped)
  y <- factor(mapped[keep], levels = c(NEGATIVE, POSITIVE))
  x <- t(expr[, keep, drop = FALSE])
  rownames(x) <- colnames(expr)[keep]
  ord <- order(rownames(x))
  list(
    x = x[ord, , drop = FALSE], y = y[ord], label_col = label_col,
    total_n = length(mapped), excluded_n = sum(!keep),
    excluded_zero = sum(bare == "0", na.rm = TRUE),
    excluded_missing_or_other = sum(!keep & bare != "0", na.rm = TRUE)
  )
}

select_discovery_cost <- function(x, y) {
  set.seed(SEED)
  folds <- caret::createFolds(y, k = CV_FOLDS, returnTrain = FALSE)
  scores <- vapply(COST_GRID, function(cost) {
    fold_scores <- numeric(length(folds))
    for (f in seq_along(folds)) {
      validation <- folds[[f]]
      training <- setdiff(seq_along(y), validation)
      features <- select_features(
        x[training, , drop = FALSE], y[training], method = "t_test", top_k = TOP_K
      )$features
      scaler <- fit_standardizer(x[training, features, drop = FALSE])
      x_tr <- apply_standardizer(x[training, features, drop = FALSE], scaler)
      x_va <- apply_standardizer(x[validation, features, drop = FALSE], scaler)
      model <- train_smo_svm(
        x_tr, y[training], kernel = KERNEL, cost = cost,
        class_weights = TRUE, probability = TRUE
      )
      probability <- predict_smo_svm(
        model, x_va, positive_class = POSITIVE
      )$probability_positive
      fold_scores[[f]] <- compute_operating_metrics(
        y[validation], probability, POSITIVE, THRESHOLD
      )$balanced_accuracy
    }
    mean(fold_scores, na.rm = TRUE)
  }, numeric(1))
  list(
    best_cost = COST_GRID[[which.max(scores)]],
    scores = data.frame(cost = COST_GRID, mean_cv_balanced_accuracy = scores)
  )
}

fit_frozen_discovery_model <- function(x, y) {
  # This ordering mirrors the original GSE25065 script: global discovery
  # feature set, discovery-only guarded cost selection, then one final fit.
  selection <- select_features(x, y, method = "t_test", top_k = TOP_K)
  features <- selection$features
  cost_selection <- select_discovery_cost(x, y)
  scaler <- fit_standardizer(x[, features, drop = FALSE])
  model <- train_smo_svm(
    apply_standardizer(x[, features, drop = FALSE], scaler), y,
    kernel = KERNEL, cost = cost_selection$best_cost,
    class_weights = TRUE, probability = TRUE
  )
  list(
    features = features,
    feature_scores = selection$scores,
    cost = cost_selection$best_cost,
    cost_scores = cost_selection$scores,
    scaler = scaler,
    model = model
  )
}

apply_frozen <- function(frozen, cohort, cohort_name, scaling, analysis_role) {
  missing <- setdiff(frozen$features, colnames(cohort$x))
  if (length(missing)) {
    stop(
      cohort_name, " is missing ", length(missing),
      " frozen features; unified exact-model projection cannot proceed."
    )
  }
  scaler <- if (identical(scaling, "discovery_derived_primary")) {
    frozen$scaler
  } else if (identical(scaling, "within_cohort_zscore_sensitivity")) {
    fit_standardizer(cohort$x[, frozen$features, drop = FALSE])
  } else {
    stop("Unknown scaling variant.")
  }
  x_scaled <- apply_standardizer(cohort$x[, frozen$features, drop = FALSE], scaler)
  prediction <- predict_smo_svm(
    frozen$model, x_scaled, positive_class = POSITIVE
  )
  probability <- prediction$probability_positive
  explicit_class <- ifelse(probability >= THRESHOLD, POSITIVE, NEGATIVE)
  data.frame(
    cohort = cohort_name,
    analysis_role = analysis_role,
    scaling = scaling,
    sample_id = rownames(cohort$x),
    truth = as.character(cohort$y),
    probability_pcr = probability,
    predicted_threshold_0_5 = explicit_class,
    predicted_package_class = prediction$predicted_class,
    package_threshold_discordant = explicit_class != prediction$predicted_class,
    stringsAsFactors = FALSE
  )
}

point_estimates <- function(predictions) {
  groups <- split(predictions, interaction(predictions$cohort, predictions$scaling,
                                           drop = TRUE))
  out <- do.call(rbind, lapply(groups, function(d) {
    metrics <- evaluate_binary_predictions(
      d$truth, d$probability_pcr, POSITIVE, THRESHOLD
    )
    data.frame(
      cohort = d$cohort[[1]], analysis_role = d$analysis_role[[1]],
      scaling = d$scaling[[1]], metrics,
      package_threshold_discordance_n = sum(d$package_threshold_discordant),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

bootstrap_metrics <- c(
  "auroc", "pr_auc", "balanced_accuracy", "mcc", "sensitivity",
  "specificity", "brier_score", "brier_skill_score", "log_loss",
  "calibration_intercept", "calibration_slope", "calibration_in_the_large"
)

conditional_bootstrap <- function(predictions, B = BOOTSTRAP_B) {
  groups <- split(predictions, interaction(predictions$cohort, predictions$scaling,
                                           drop = TRUE))
  rows <- list()
  group_index <- 0L
  for (d in groups) {
    group_index <- group_index + 1L
    pos <- which(d$truth == POSITIVE)
    neg <- which(d$truth == NEGATIVE)
    point <- evaluate_binary_predictions(
      d$truth, d$probability_pcr, POSITIVE, THRESHOLD
    )
    boot <- matrix(NA_real_, nrow = B, ncol = length(bootstrap_metrics),
                   dimnames = list(NULL, bootstrap_metrics))
    set.seed(SEED + 600000L + group_index)
    for (b in seq_len(B)) {
      idx <- c(sample(pos, replace = TRUE), sample(neg, replace = TRUE))
      bm <- evaluate_binary_predictions(
        d$truth[idx], d$probability_pcr[idx], POSITIVE, THRESHOLD
      )
      boot[b, ] <- as.numeric(bm[1, bootstrap_metrics])
    }
    for (metric in bootstrap_metrics) {
      ci <- stats::quantile(boot[, metric], c(0.025, 0.975),
                            na.rm = TRUE, names = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = d$cohort[[1]], analysis_role = d$analysis_role[[1]],
        scaling = d$scaling[[1]], metric = metric,
        point = point[[metric]], ci_lo = ci[[1]], ci_hi = ci[[2]],
        n_boot = B,
        method = paste(
          "stratified percentile bootstrap conditional on the single frozen",
          "model and stored predictions; discovery workflow not refit"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

write_external_figure <- function(points, intervals) {
  variants <- data.frame(
    cohort = c("GSE25065", "GSE41998", "GSE41998"),
    scaling = c("discovery_derived_primary", "discovery_derived_primary",
                "within_cohort_zscore_sensitivity"),
    label = c("GSE25065\nprimary", "GSE41998\nprimary", "GSE41998\nz-score"),
    stringsAsFactors = FALSE
  )
  pick_points <- do.call(rbind, lapply(seq_len(nrow(variants)), function(i) {
    points[points$cohort == variants$cohort[[i]] &
             points$scaling == variants$scaling[[i]], ]
  }))
  pick_interval <- function(i, metric) {
    intervals[intervals$cohort == variants$cohort[[i]] &
                intervals$scaling == variants$scaling[[i]] &
                intervals$metric == metric, ]
  }
  for (device in list(
    function() grDevices::png(file.path(FIGURES_DIR, "external_validation_unified.png"),
                              width = 1650, height = 620, res = 120),
    function() grDevices::pdf(file.path(FIGURES_DIR, "external_validation_unified.pdf"),
                              width = 13.7, height = 5.1)
  )) {
    device()
    graphics::par(mfrow = c(1, 3), mar = c(5.2, 4.2, 3, 1))

    x <- seq_len(nrow(variants))
    graphics::plot(x, pick_points$auroc, pch = 19, cex = 1.15, xaxt = "n",
                   ylim = c(0.2, 0.85), xlab = "", ylab = "Metric value",
                   main = "Discrimination")
    for (i in x) {
      ci <- pick_interval(i, "auroc")
      graphics::arrows(i, ci$ci_lo, i, ci$ci_hi, angle = 90, code = 3,
                       length = 0.05, col = "#1D4ED8")
    }
    graphics::points(x, pick_points$pr_auc, pch = 17, cex = 1.1, col = "#7C3AED")
    for (i in x) {
      ci <- pick_interval(i, "pr_auc")
      graphics::arrows(i + 0.05, ci$ci_lo, i + 0.05, ci$ci_hi, angle = 90,
                       code = 3, length = 0.05, col = "#7C3AED")
    }
    graphics::axis(1, at = x, labels = variants$label, cex.axis = 0.82)
    graphics::legend("bottomright", bty = "n", pch = c(19, 17),
                     col = c("#1D4ED8", "#7C3AED"), legend = c("AUROC", "PR-AUC"))

    op <- rbind(sensitivity = pick_points$sensitivity,
                specificity = pick_points$specificity)
    graphics::barplot(op, beside = TRUE, ylim = c(0, 1),
                      names.arg = variants$label, cex.names = 0.82,
                      col = c("#F59E0B", "#0F766E"), border = NA,
                      ylab = "Proportion", main = "Operating point at P(pCR) = 0.50")
    graphics::legend("top", bty = "n", fill = c("#F59E0B", "#0F766E"),
                     legend = c("Sensitivity", "Specificity"), horiz = TRUE,
                     cex = 0.8)

    slope <- pick_points$calibration_slope
    graphics::plot(x, slope, pch = 19, cex = 1.15, xaxt = "n", ylim = c(-0.2, 1.2),
                   xlab = "", ylab = "Calibration slope",
                   main = "Calibration transport")
    for (i in x) {
      ci <- pick_interval(i, "calibration_slope")
      graphics::arrows(i, ci$ci_lo, i, ci$ci_hi, angle = 90, code = 3,
                       length = 0.05, col = "#B91C1C")
    }
    graphics::abline(h = 1, lty = 3)
    graphics::axis(1, at = x, labels = variants$label, cex.axis = 0.82)
    graphics::mtext("Bars: conditional 95% bootstrap intervals", side = 1,
                    line = 4.1, cex = 0.75)
    grDevices::dev.off()
  }
}

legacy_regression_checks <- function(points) {
  anchors <- data.frame(
    cohort = c("GSE25065", "GSE25065", "GSE41998", "GSE41998", "GSE41998"),
    scaling = c(
      "discovery_derived_primary", "discovery_derived_primary",
      "discovery_derived_primary", "discovery_derived_primary",
      "within_cohort_zscore_sensitivity"
    ),
    metric = c("auroc", "pr_auc", "auroc", "pr_auc", "auroc"),
    legacy_value = c(0.5633, 0.2986, 0.6922, 0.4634, 0.7067),
    tolerance = 0.005,
    stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(seq_len(nrow(anchors)), function(i) {
    a <- anchors[i, ]
    row <- points[points$cohort == a$cohort & points$scaling == a$scaling, ]
    current <- row[[a$metric]]
    data.frame(
      a,
      unified_value = current,
      absolute_difference = abs(current - a$legacy_value),
      pass = abs(current - a$legacy_value) <= a$tolerance,
      stringsAsFactors = FALSE
    )
  }))
}

write_cohort_context <- function(gse25065, gse41998) {
  context <- data.frame(
    cohort = c("GSE25055", "GSE25065", "GSE41998"),
    role = c(
      "discovery/model development", "same-study-family external validation",
      "independent cross-platform transportability sensitivity"
    ),
    platform = c("GPL96 / HG-U133A", "GPL96 / HG-U133A", "GPL571 / HG-U133A 2.0"),
    analyzed_n = c(306L, length(gse25065$y), length(gse41998$y)),
    pcr_n = c(57L, sum(gse25065$y == POSITIVE), sum(gse41998$y == POSITIVE)),
    rd_n = c(249L, sum(gse25065$y == NEGATIVE), sum(gse41998$y == NEGATIVE)),
    excluded_n = c(4L, gse25065$excluded_n, gse41998$excluded_n),
    endpoint_mapping = c(
      "pathologic_response_pcr_rd: pCR/RD; NA excluded",
      "pathologic_response_pcr_rd: pCR/RD; NA excluded",
      "pcr: Yes->pCR, No->RD; literal 0 and missing/other excluded"
    ),
    primary_scaling = c(
      "training-partition/discovery-derived", "frozen discovery-derived",
      "frozen discovery-derived"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(context, file.path(TABLES_DIR, "cohort_context.csv"), row.names = FALSE)
  context
}

main <- function() {
  start <- Sys.time()
  if (!requireNamespace("caret", quietly = TRUE)) stop("caret required.")
  discovery <- load_discovery()
  message(sprintf("[unified] fitting one frozen model: discovery n=%d p=%d",
                  nrow(discovery$x), ncol(discovery$x)))
  frozen <- fit_frozen_discovery_model(discovery$x, discovery$y)
  message(sprintf("[unified] frozen top-%d, cost=%s; now loading external cohorts.",
                  length(frozen$features), frozen$cost))

  # External cohorts are loaded only after the discovery pipeline is frozen.
  gse25065 <- load_binary_gse("GSE25065")
  gse41998 <- load_gse41998()
  assert_no_overlap(rownames(discovery$x), rownames(gse25065$x))
  assert_no_overlap(rownames(discovery$x), rownames(gse41998$x))

  predictions <- rbind(
    apply_frozen(
      frozen, gse25065, "GSE25065", "discovery_derived_primary",
      "same-study-family external validation"
    ),
    apply_frozen(
      frozen, gse41998, "GSE41998", "discovery_derived_primary",
      "cross-platform transportability sensitivity"
    ),
    apply_frozen(
      frozen, gse41998, "GSE41998", "within_cohort_zscore_sensitivity",
      "label-blind scaling sensitivity"
    )
  )
  points <- point_estimates(predictions)
  intervals <- conditional_bootstrap(predictions)
  write_external_figure(points, intervals)
  checks <- legacy_regression_checks(points)
  context <- write_cohort_context(gse25065, gse41998)

  feature_manifest <- data.frame(
    rank = seq_along(frozen$features), feature = frozen$features,
    abs_t = as.numeric(frozen$feature_scores),
    present_gse25065 = frozen$features %in% colnames(gse25065$x),
    present_gse41998 = frozen$features %in% colnames(gse41998$x),
    stringsAsFactors = FALSE
  )
  model_manifest <- data.frame(
    discovery = DISCOVERY,
    discovery_n = nrow(discovery$x),
    top_k = TOP_K,
    svm_kernel = KERNEL,
    selected_cost = frozen$cost,
    class_weights = "inverse class frequency from discovery",
    probability = "e1071/libsvm probability=TRUE; one final fit reused",
    primary_scaling = "GSE25055 center/scale frozen and reused",
    threshold = THRESHOLD,
    stringsAsFactors = FALSE
  )

  write.csv(predictions, file.path(RESULTS_DIR, "external_predictions.csv"), row.names = FALSE)
  write.csv(feature_manifest, file.path(RESULTS_DIR, "frozen_feature_manifest.csv"), row.names = FALSE)
  write.csv(frozen$cost_scores, file.path(RESULTS_DIR, "discovery_cost_selection.csv"), row.names = FALSE)
  write.csv(model_manifest, file.path(RESULTS_DIR, "frozen_model_manifest.csv"), row.names = FALSE)
  write.csv(points, file.path(TABLES_DIR, "external_point_estimates.csv"), row.names = FALSE)
  write.csv(intervals, file.path(TABLES_DIR, "external_bootstrap_intervals.csv"), row.names = FALSE)
  write.csv(checks, file.path(RESULTS_DIR, "legacy_anchor_checks.csv"), row.names = FALSE)

  p65 <- points[points$cohort == "GSE25065", ]
  p98 <- points[points$cohort == "GSE41998" &
                  points$scaling == "discovery_derived_primary", ]
  s98 <- points[points$cohort == "GSE41998" &
                  points$scaling == "within_cohort_zscore_sensitivity", ]
  notes <- c(
    "# Unified frozen-model external validation", "",
    sprintf("- One GSE25055 model was fit once and reused for all projections (top-%d; cost=%s).",
            TOP_K, frozen$cost),
    sprintf("- Both primary external projections use the identical discovery-derived scaler and explicit P(pCR) >= %.2f rule.", THRESHOLD),
    sprintf("- Package-class versus explicit-threshold discordance: %d predictions.",
            sum(predictions$package_threshold_discordant)),
    sprintf("- Frozen feature recovery: GSE25065 %d/%d; GSE41998 %d/%d.",
            sum(feature_manifest$present_gse25065), TOP_K,
            sum(feature_manifest$present_gse41998), TOP_K),
    "", "## Primary external results",
    sprintf("- GSE25065 (n=%d; pCR prevalence %.3f): AUROC %.4f, PR-AUC %.4f, Brier %.4f, calibration intercept %.3f, slope %.3f; sensitivity %.3f, specificity %.3f.",
            p65$n, p65$prevalence, p65$auroc, p65$pr_auc, p65$brier_score,
            p65$calibration_intercept, p65$calibration_slope,
            p65$sensitivity, p65$specificity),
    sprintf("- GSE41998 primary (n=%d; pCR prevalence %.3f): AUROC %.4f, PR-AUC %.4f, Brier %.4f, calibration intercept %.3f, slope %.3f; sensitivity %.3f, specificity %.3f.",
            p98$n, p98$prevalence, p98$auroc, p98$pr_auc, p98$brier_score,
            p98$calibration_intercept, p98$calibration_slope,
            p98$sensitivity, p98$specificity),
    sprintf("- GSE41998 label-blind z-score sensitivity: AUROC %.4f, PR-AUC %.4f, Brier %.4f; sensitivity %.3f, specificity %.3f.",
            s98$auroc, s98$pr_auc, s98$brier_score,
            s98$sensitivity, s98$specificity),
    "", "## Interpretation and uncertainty",
    "Threshold-independent discrimination is primary. PR-AUC must be interpreted against each cohort's pCR prevalence.",
    "The same 0.5 operating point behaves very differently across cohorts and scaling variants; this is evidence of transport/calibration instability, not a reason to retune on external outcomes.",
    "Bootstrap intervals are conditional on the single frozen model and stored external predictions. They quantify external sample variation, not discovery-model refitting uncertainty.",
    "GSE41998 combines platform, treatment, population, and endpoint-context transport; its result must not be attributed to platform alone.",
    sprintf("Legacy rounded-result regression checks: %d/%d passed within tolerance.",
            sum(checks$pass), nrow(checks)),
    "", sprintf("_Runtime %.1f minutes._", as.numeric(difftime(Sys.time(), start, units = "mins")))
  )
  writeLines(notes, file.path(RESULTS_DIR, "external_validation_unified_notes.md"))
  if (!all(checks$pass)) warning("One or more legacy rounded-result checks failed; inspect before claim updates.")
  message(sprintf("[unified] complete in %.1f minutes; legacy checks %d/%d.",
                  as.numeric(difftime(Sys.time(), start, units = "mins")),
                  sum(checks$pass), nrow(checks)))
  invisible(list(points = points, intervals = intervals, checks = checks, context = context))
}

if (sys.nframe() == 0L) main()
