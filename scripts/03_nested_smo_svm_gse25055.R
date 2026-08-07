# =============================================================================
# scripts/03_nested_smo_svm_gse25055.R
# -----------------------------------------------------------------------------
# Pilot Pipeline B: GUARDED NESTED SMO/SVM pipeline on GSE25055 ONLY.
#
# This is the leakage-controlled counterpart to script 02. Every supervised,
# data-dependent step happens INSIDE training folds only:
#
#   * Feature selection (t-test, top K = 100) -> fit on inner/outer TRAIN only.
#   * Standardization (center/scale)          -> fit on TRAIN only, applied to
#                                                 validation/test folds.
#   * Hyperparameter tuning (SVM cost)         -> chosen by INNER CV only.
#   * Class imbalance                          -> handled via class WEIGHTS
#                                                 (no SMOTE in the first pilot).
#
# Outer CV gives the honest generalization estimate; inner CV does selection
# and tuning. The outer test fold is NEVER seen during any fitting step.
#
# Scope (pilot):
#   - GSE25055 ONLY. Do NOT load GSE25065 / GSE41998 / GSE20194 / GSE20271.
#   - Labels: pCR vs RD from 'pathologic_response_pcr_rd'; NA-coded excluded.
# =============================================================================

# ---- Setup ------------------------------------------------------------------
source("R/00_config.R")          # SEED, %||%
source("R/feature_selection.R")  # select_features()
source("R/preprocessing.R")      # filter_near_zero_variance(), standardizer
source("R/model_smo_svm.R")      # train_smo_svm(), predict_smo_svm()
source("R/metrics.R")            # compute_binary_metrics()
source("R/leakage_checks.R")     # assert_no_overlap(), log_pipeline_step()

set.seed(SEED)

# Pilot configuration -----------------------------------------------------
ACCESSION   <- "GSE25055"
LABEL_FIELD <- "pathologic_response_pcr_rd"
POSITIVE    <- "pCR"
TOP_K       <- 100                 # first pilot feature budget
OUTER_FOLDS <- 5
INNER_FOLDS <- 5
COST_GRID   <- c(0.25, 1, 4)       # small linear-SVM cost grid for the pilot
KERNEL      <- "linear"

RESULTS_DIR <- file.path("results", "pilot_gse25055")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Data loading (identical contract to script 02) -------------------------
# Reuses the same loader logic so both pipelines see EXACTLY the same samples,
# labels, and NA exclusions. This keeps the leaky-vs-guarded comparison fair.
load_gse25055 <- function(accession = ACCESSION, label_field = LABEL_FIELD) {
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
    stop("Unexpected label values after NA exclusion: ",
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

# Label-independent: removing near-zero-variance probes uses no labels, so it
# may be applied once up front without leakage.
prep_expression <- function(x) filter_near_zero_variance(x, cutoff = 1e-8)$x

# ---- Inner CV: pick the best SVM cost on TRAINING data only ------------------
# GUARDED: feature selection + scaling + tuning all happen on inner-train folds
# carved from the OUTER-TRAIN data. The outer-test fold is untouched here.
tune_inner <- function(x_train, y_train, top_k = TOP_K, cost_grid = COST_GRID,
                       inner_folds = INNER_FOLDS, positive = POSITIVE) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package caret is required for fold creation.")
  }
  set.seed(SEED)  # deterministic inner-fold assignment
  folds <- caret::createFolds(y_train, k = inner_folds, returnTrain = FALSE)

  score_for_cost <- function(cost) {
    bal_accs <- numeric(0)
    for (f in seq_along(folds)) {
      val_idx <- folds[[f]]
      tr_idx  <- setdiff(seq_along(y_train), val_idx)

      # ---- LEAKAGE GUARD: select features on inner-TRAIN only ----------------
      fs   <- select_features(x_train[tr_idx, , drop = FALSE], y_train[tr_idx],
                              method = "t_test", top_k = top_k)
      feats <- fs$features

      # ---- LEAKAGE GUARD: fit scaler on inner-TRAIN only ---------------------
      std  <- fit_standardizer(x_train[tr_idx, feats, drop = FALSE])
      x_tr <- apply_standardizer(x_train[tr_idx,  feats, drop = FALSE], std)
      x_va <- apply_standardizer(x_train[val_idx, feats, drop = FALSE], std)

      model <- train_smo_svm(x_tr, y_train[tr_idx], kernel = KERNEL,
                             cost = cost, class_weights = TRUE,
                             probability = TRUE)
      pr <- predict_smo_svm(model, x_va, positive_class = positive)
      sens <- sum(y_train[val_idx] == positive & pr$predicted_class == positive) /
              sum(y_train[val_idx] == positive)
      spec <- sum(y_train[val_idx] != positive & pr$predicted_class != positive) /
              sum(y_train[val_idx] != positive)
      bal_accs <- c(bal_accs, mean(c(sens, spec), na.rm = TRUE))
    }
    mean(bal_accs, na.rm = TRUE)
  }

  scores    <- vapply(cost_grid, score_for_cost, numeric(1))
  best_cost <- cost_grid[which.max(scores)]
  log_pipeline_step("inner_tune", notes = sprintf("best_cost=%s (bal_acc=%.3f)",
                                                   best_cost, max(scores)))
  best_cost
}

# ---- Outer CV: honest performance estimate ----------------------------------
run_nested <- function(x, y, outer_folds = OUTER_FOLDS, top_k = TOP_K,
                       positive = POSITIVE) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package caret is required for fold creation.")
  }
  set.seed(SEED)
  folds <- caret::createFolds(y, k = outer_folds, returnTrain = FALSE)

  preds        <- list()
  fold_records <- list()
  feat_records <- list()

  for (f in seq_along(folds)) {
    test_idx  <- folds[[f]]
    train_idx <- setdiff(seq_along(y), test_idx)

    # Defensive leakage check: outer train/test must not share samples.
    assert_no_overlap(rownames(x)[train_idx], rownames(x)[test_idx])

    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]

    # ---- INNER CV: tune cost using outer-TRAIN only ----------------------
    best_cost <- tune_inner(x_train, y_train, top_k = top_k)

    # ---- Refit on full outer-TRAIN with guarded steps -------------------
    # Feature selection on outer-TRAIN only (NOT on the full dataset).
    fs    <- select_features(x_train, y_train, method = "t_test", top_k = top_k)
    feats <- fs$features

    # Record per-outer-fold selected features for feature-stability analysis.
    # Pure logging: does not alter selection, model logic, top K, folds, or tuning.
    feat_records[[f]] <- data.frame(
      fold    = f,
      rank    = seq_along(feats),
      feature = feats,
      stringsAsFactors = FALSE
    )

    # Scaler fit on outer-TRAIN only, then applied to the outer-TEST fold.
    std   <- fit_standardizer(x_train[, feats, drop = FALSE])
    x_tr  <- apply_standardizer(x_train[, feats, drop = FALSE], std)
    x_te  <- apply_standardizer(x[test_idx, feats, drop = FALSE], std)

    model <- train_smo_svm(x_tr, y_train, kernel = KERNEL, cost = best_cost,
                           class_weights = TRUE, probability = TRUE)
    pr <- predict_smo_svm(model, x_te, positive_class = positive)

    preds[[f]] <- data.frame(
      fold      = f,
      sample_id = rownames(x)[test_idx],
      truth     = as.character(y[test_idx]),
      predicted = pr$predicted_class,
      prob_pos  = pr$probability_positive,
      stringsAsFactors = FALSE
    )
    fold_records[[f]] <- data.frame(
      fold        = f,
      best_cost   = best_cost,
      n_features  = length(feats),
      stringsAsFactors = FALSE
    )
  }

  list(
    predictions = do.call(rbind, preds),
    fold_info   = do.call(rbind, fold_records),
    selected_features = do.call(rbind, feat_records)
  )
}

# ---- Summarize pooled outer-fold predictions into metrics -------------------
summarize_metrics <- function(pred_df, positive = POSITIVE) {
  m <- compute_binary_metrics(truth = pred_df$truth, estimate = pred_df$predicted,
                              probability = pred_df$prob_pos, positive = positive)
  sens <- sum(pred_df$truth == positive & pred_df$predicted == positive) /
          sum(pred_df$truth == positive)
  spec <- sum(pred_df$truth != positive & pred_df$predicted != positive) /
          sum(pred_df$truth != positive)
  data.frame(
    pipeline          = "B_guarded_nested",
    auroc             = m$auroc,
    pr_auc            = m$pr_auc,
    balanced_accuracy = m$balanced_accuracy,
    mcc               = m$mcc,
    sensitivity       = sens,
    specificity       = spec,
    stringsAsFactors  = FALSE
  )
}

# ---- Main -------------------------------------------------------------------
main <- function() {
  dat <- load_gse25055()
  x   <- prep_expression(dat$x)
  y   <- dat$y

  message(sprintf("[main] Class balance: %s",
                  paste(names(table(y)), table(y), sep = "=", collapse = ", ")))

  res     <- run_nested(x, y, top_k = TOP_K)
  metrics <- summarize_metrics(res$predictions)

  utils::write.csv(metrics,
                   file.path(RESULTS_DIR, "nested_smo_svm_metrics.csv"),
                   row.names = FALSE)
  utils::write.csv(res$predictions,
                   file.path(RESULTS_DIR, "nested_smo_svm_predictions.csv"),
                   row.names = FALSE)
  utils::write.csv(res$fold_info,
                   file.path(RESULTS_DIR, "nested_smo_svm_fold_info.csv"),
                   row.names = FALSE)
  utils::write.csv(res$selected_features,
                   file.path(RESULTS_DIR, "nested_selected_features_by_fold.csv"),
                   row.names = FALSE)

  message("[main] Guarded nested pipeline complete. Metrics:")
  print(metrics)
  invisible(list(metrics = metrics, predictions = res$predictions,
                 fold_info = res$fold_info))
}

if (sys.nframe() == 0) {
  main()
}
