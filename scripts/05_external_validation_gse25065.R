# =============================================================================
# scripts/05_external_validation_gse25065.R
# -----------------------------------------------------------------------------
# Task 5: First EXTERNAL VALIDATION analysis.
#
#   Discovery / training cohort : GSE25055  (Hatzis et al., JAMA 2011)
#   External validation cohort  : GSE25065  (Hatzis et al., JAMA 2011)
#
# GSE25065 is treated as a TRUE external validation cohort. It is NOT used for:
#   * feature selection
#   * hyperparameter (SVM cost) tuning
#   * scaling-parameter (center/scale) estimation
#   * decision-threshold selection
#   * class-balancing decisions
#   * any model-design decision
#
# The full pipeline (feature set, scaling parameters, SVM model, decision rule)
# is FIT and FROZEN on GSE25055 only, then applied exactly once to GSE25065.
#
# Constraints honored: GSE25055 + GSE25065 only (the two are non-overlapping by
# design per docs/dataset_audit_report.md). No GSE41998 / GSE20194 / GSE20271.
# No raw CEL files. Series-matrix metadata only, via GEOquery. Same t-test
# top K = 100 strategy, same linear kernel, and same cost grid as the pilot.
# =============================================================================

# ---- Setup ------------------------------------------------------------------
source("R/00_config.R")          # SEED, %||%
source("R/feature_selection.R")  # select_features()
source("R/preprocessing.R")      # filter_near_zero_variance(), standardizer, align_features()
source("R/model_smo_svm.R")      # train_smo_svm(), predict_smo_svm()
source("R/metrics.R")            # compute_binary_metrics()
source("R/leakage_checks.R")     # assert_no_overlap(), assert_features_available(), log_pipeline_step()

set.seed(SEED)

# ---- Configuration ----------------------------------------------------------
DISCOVERY   <- "GSE25055"
EXTERNAL    <- "GSE25065"
LABEL_FIELD <- "pathologic_response_pcr_rd"
POSITIVE    <- "pCR"
TOP_K       <- 100                 # identical feature budget to the pilot
CV_FOLDS    <- 5                   # discovery-only CV for SVM-cost selection
COST_GRID   <- c(0.25, 1, 4)       # identical cost grid to the pilot
KERNEL      <- "linear"

RESULTS_DIR <- file.path("results", "external_validation_gse25065")
TABLES_DIR  <- file.path("tables",  "external_validation_gse25065")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)

# ---- Generic GEO loader (same contract as the pilot's load_gse25055) --------
# Parameterized by accession so discovery and external cohorts are loaded by
# IDENTICAL logic: same label field, same NA handling, same factor levels.
load_gse <- function(accession, label_field = LABEL_FIELD) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package GEOquery (Bioconductor) is required to load ", accession, ".")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package Biobase (Bioconductor) is required.")
  }
  gset  <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)
  eset  <- gset[[1]]
  expr  <- Biobase::exprs(eset)
  pheno <- Biobase::pData(eset)

  label_col <- NULL
  for (cn in colnames(pheno)) {
    vals <- as.character(pheno[[cn]])
    if (any(grepl(label_field, vals, ignore.case = TRUE)) ||
        grepl(label_field, cn, ignore.case = TRUE)) {
      label_col <- cn; break
    }
  }
  if (is.null(label_col)) {
    stop("Could not find label field '", label_field, "' in ", accession, ".")
  }

  raw_labels <- as.character(pheno[[label_col]])
  raw_labels <- sub(paste0("^.*", label_field, "\\s*[:=]?\\s*"), "", raw_labels,
                    ignore.case = TRUE)
  raw_labels <- trimws(raw_labels)

  is_na <- is.na(raw_labels) | raw_labels %in% c("NA", "na", "N/A", "", "NaN")
  keep  <- !is_na
  message(sprintf("[load] %s: %d total; excluding %d NA; keeping %d.",
                  accession, length(raw_labels), sum(is_na), sum(keep)))

  labels <- factor(raw_labels[keep], levels = c("RD", "pCR"))
  if (any(is.na(labels))) {
    stop("Unexpected label values after NA exclusion in ", accession, ": ",
         paste(unique(raw_labels[keep]), collapse = ", "))
  }

  x <- t(expr[, keep, drop = FALSE])
  ids <- colnames(expr)[keep]
  rownames(x) <- ids
  ord <- order(ids)
  x <- x[ord, , drop = FALSE]
  labels <- labels[ord]
  ids <- ids[ord]
  list(x = x, y = labels, ids = ids)
}

# ---- Discovery-only SVM-cost selection (guarded CV) -------------------------
# Pure discovery-side procedure: GSE25065 is never touched. Within each CV fold
# carved from GSE25055, feature selection AND scaling are fit on the fold-train
# split only, so the cost choice is not optimistically biased.
select_cost_discovery <- function(x, y, top_k = TOP_K, cost_grid = COST_GRID,
                                  cv_folds = CV_FOLDS, positive = POSITIVE) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package caret is required for fold creation.")
  }
  set.seed(SEED)
  folds <- caret::createFolds(y, k = cv_folds, returnTrain = FALSE)

  score_for_cost <- function(cost) {
    bal <- numeric(0)
    for (f in seq_along(folds)) {
      va <- folds[[f]]
      tr <- setdiff(seq_along(y), va)
      fs    <- select_features(x[tr, , drop = FALSE], y[tr],
                               method = "t_test", top_k = top_k)
      feats <- fs$features
      std   <- fit_standardizer(x[tr, feats, drop = FALSE])
      x_tr  <- apply_standardizer(x[tr, feats, drop = FALSE], std)
      x_va  <- apply_standardizer(x[va, feats, drop = FALSE], std)
      model <- train_smo_svm(x_tr, y[tr], kernel = KERNEL, cost = cost,
                             class_weights = TRUE, probability = TRUE)
      pr    <- predict_smo_svm(model, x_va, positive_class = positive)
      sens  <- sum(y[va] == positive & pr$predicted_class == positive) /
               sum(y[va] == positive)
      spec  <- sum(y[va] != positive & pr$predicted_class != positive) /
               sum(y[va] != positive)
      bal   <- c(bal, mean(c(sens, spec), na.rm = TRUE))
    }
    mean(bal, na.rm = TRUE)
  }

  scores    <- vapply(cost_grid, score_for_cost, numeric(1))
  best_cost <- cost_grid[which.max(scores)]
  log_pipeline_step("discovery_cost_select",
                    notes = sprintf("best_cost=%s (bal_acc=%.3f)",
                                    best_cost, max(scores)))
  list(best_cost = best_cost, cost_grid = cost_grid, cv_bal_acc = scores)
}

# ---- Build the FROZEN discovery pipeline -------------------------------------
fit_frozen_pipeline <- function(x, y, top_k = TOP_K, positive = POSITIVE) {
  # 1. Feature selection on the FULL discovery cohort (no external data).
  fs    <- select_features(x, y, method = "t_test", top_k = top_k)
  feats <- fs$features

  # 2. Cost selection by discovery-only guarded CV.
  cost_sel <- select_cost_discovery(x, y, top_k = top_k, positive = positive)

  # 3. Scaling parameters fit on full discovery (selected features only).
  std <- fit_standardizer(x[, feats, drop = FALSE])

  # 4. Final SVM/SMO trained on full discovery; decision rule = libsvm default
  #    (0.5 on calibrated P(pCR)); NO post-hoc threshold tuning. Everything
  #    below is now FROZEN before GSE25065 is loaded for prediction.
  x_disc <- apply_standardizer(x[, feats, drop = FALSE], std)
  model  <- train_smo_svm(x_disc, y, kernel = KERNEL, cost = cost_sel$best_cost,
                          class_weights = TRUE, probability = TRUE)

  list(features = feats, scores = fs$scores, standardizer = std,
       model = model, best_cost = cost_sel$best_cost,
       cost_select = cost_sel, positive = positive)
}

# ---- Apply frozen pipeline to the external cohort (ONCE) --------------------
apply_to_external <- function(frozen, x_ext, y_ext, positive = POSITIVE) {
  feats <- frozen$features

  # Feature/probe alignment. GSE25055 and GSE25065 share platform GPL96, so
  # the frozen probe set is expected to be fully present in GSE25065.
  present <- intersect(feats, colnames(x_ext))
  missing <- setdiff(feats, colnames(x_ext))
  log_pipeline_step("external_align",
                    notes = sprintf("%d/%d frozen features present; %d missing.",
                                    length(present), length(feats), length(missing)))
  if (length(missing) > 0) {
    stop("Frozen features absent from ", EXTERNAL, " (cannot apply frozen model ",
         "without refitting): ", paste(utils::head(missing, 10), collapse = ", "),
         if (length(missing) > 10) " ..." else "")
  }

  # Apply FROZEN scaling (GSE25055-derived center/scale). No refit on GSE25065.
  x_ext_std <- apply_standardizer(x_ext[, feats, drop = FALSE], frozen$standardizer)

  pr <- predict_smo_svm(frozen$model, x_ext_std, positive_class = positive)

  predictions <- data.frame(
    sample_id = rownames(x_ext),
    truth     = as.character(y_ext),
    predicted = pr$predicted_class,
    prob_pos  = pr$probability_positive,
    stringsAsFactors = FALSE
  )
  list(predictions = predictions, n_features_aligned = length(present),
       n_features_missing = length(missing))
}

# ---- Metrics + confusion matrix on the external cohort ----------------------
external_metrics <- function(pred_df, positive = POSITIVE) {
  m <- compute_binary_metrics(truth = pred_df$truth, estimate = pred_df$predicted,
                              probability = pred_df$prob_pos, positive = positive)
  tp <- sum(pred_df$truth == positive & pred_df$predicted == positive)
  fn <- sum(pred_df$truth == positive & pred_df$predicted != positive)
  tn <- sum(pred_df$truth != positive & pred_df$predicted != positive)
  fp <- sum(pred_df$truth != positive & pred_df$predicted == positive)
  sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  data.frame(
    analysis          = "external_validation",
    discovery_cohort  = DISCOVERY,
    external_cohort   = EXTERNAL,
    n_external        = nrow(pred_df),
    auroc             = m$auroc,
    pr_auc            = m$pr_auc,
    balanced_accuracy = m$balanced_accuracy,
    mcc               = m$mcc,
    sensitivity       = sens,
    specificity       = spec,
    tp = tp, fp = fp, tn = tn, fn = fn,
    stringsAsFactors  = FALSE
  )
}

# ---- Main -------------------------------------------------------------------
main <- function() {
  # --- Load discovery (GSE25055) ---
  disc   <- load_gse(DISCOVERY)
  x_disc <- filter_near_zero_variance(disc$x, cutoff = 1e-8)$x
  y_disc <- disc$y
  message(sprintf("[main] %s class balance: %s", DISCOVERY,
                  paste(names(table(y_disc)), table(y_disc), sep = "=", collapse = ", ")))

  # --- Fit & FREEZE the full pipeline on discovery ONLY ---
  frozen <- fit_frozen_pipeline(x_disc, y_disc, top_k = TOP_K)
  message(sprintf("[main] Frozen pipeline: %d features, cost=%s.",
                  length(frozen$features), frozen$best_cost))

  # --- Load external (GSE25065) AFTER the pipeline is frozen ---
  ext   <- load_gse(EXTERNAL)
  x_ext <- ext$x          # raw probe matrix; no GSE25065-derived filtering
  y_ext <- ext$y
  message(sprintf("[main] %s class balance: %s", EXTERNAL,
                  paste(names(table(y_ext)), table(y_ext), sep = "=", collapse = ", ")))

  # Defensive overlap check: discovery and external GSM ids must be disjoint.
  assert_no_overlap(rownames(x_disc), rownames(x_ext))

  # --- Apply frozen model to external cohort ONCE ---
  applied <- apply_to_external(frozen, x_ext, y_ext)
  metrics <- external_metrics(applied$predictions)

  message("[main] External-validation metrics on ", EXTERNAL, ":")
  print(metrics)

  # --- Write outputs ---
  utils::write.csv(applied$predictions,
                   file.path(RESULTS_DIR, "gse25065_external_predictions.csv"),
                   row.names = FALSE)
  utils::write.csv(metrics,
                   file.path(RESULTS_DIR, "gse25065_external_metrics.csv"),
                   row.names = FALSE)

  feat_tab <- data.frame(
    rank    = seq_along(frozen$features),
    feature = frozen$features,
    abs_t   = as.numeric(frozen$scores),
    stringsAsFactors = FALSE
  )
  utils::write.csv(feat_tab,
                   file.path(RESULTS_DIR, "final_gse25055_selected_features.csv"),
                   row.names = FALSE)

  # Compact key/value summary table.
  summary_tab <- data.frame(
    metric = c("discovery_cohort", "external_cohort",
               "discovery_n", "discovery_pCR", "discovery_RD",
               "external_n", "external_pCR", "external_RD",
               "top_k", "selected_features", "features_aligned",
               "features_missing", "svm_kernel", "svm_cost",
               "auroc", "pr_auc", "balanced_accuracy", "mcc",
               "sensitivity", "specificity",
               "tp", "fp", "tn", "fn"),
    value = c(DISCOVERY, EXTERNAL,
              nrow(x_disc), sum(y_disc == "pCR"), sum(y_disc == "RD"),
              nrow(x_ext), sum(y_ext == "pCR"), sum(y_ext == "RD"),
              TOP_K, length(frozen$features), applied$n_features_aligned,
              applied$n_features_missing, KERNEL, frozen$best_cost,
              round(metrics$auroc, 4), round(metrics$pr_auc, 4),
              round(metrics$balanced_accuracy, 4), round(metrics$mcc, 4),
              round(metrics$sensitivity, 4), round(metrics$specificity, 4),
              metrics$tp, metrics$fp, metrics$tn, metrics$fn),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary_tab,
                   file.path(TABLES_DIR, "external_validation_summary.csv"),
                   row.names = FALSE)

  message("[main] External validation complete. Outputs written to ",
          RESULTS_DIR, " and ", TABLES_DIR, ".")
  invisible(list(frozen = frozen, predictions = applied$predictions,
                 metrics = metrics, summary = summary_tab))
}

if (sys.nframe() == 0) {
  res05 <- main()
}
