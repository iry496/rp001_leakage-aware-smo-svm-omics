# =============================================================================
# scripts/02_leaky_baseline_gse25055.R
# -----------------------------------------------------------------------------
# Pilot Pipeline A: NAIVE / LEAKY SMO/SVM baseline on GSE25055 ONLY.
#
# >>> WARNING <<<
# This script is INTENTIONALLY LEAKY. It performs supervised feature selection
# on the ENTIRE dataset BEFORE cross-validation. This deliberately leaks label
# information from the held-out folds into feature selection and is included
# ONLY as a comparison ("how optimistic does naive practice look?").
# DO NOT reuse this pattern in the guarded pipeline (script 03) or anywhere
# results are reported as honest generalization estimates.
#
# Scope (pilot):
#   - GSE25055 ONLY. Do NOT load GSE25065 / GSE41998 / GSE20194 / GSE20271.
#   - Labels: pCR vs RD from the 'pathologic_response_pcr_rd' metadata field.
#   - Samples coded NA are excluded.
#   - Simple t-test feature selector, top K = 100.
#   - Linear SVM (libsvm via e1071), performance via cross-validation.
# =============================================================================

# ---- Setup ------------------------------------------------------------------
source("R/00_config.R")          # defines SEED, %||%, etc.
source("R/feature_selection.R")  # select_features() / select_ttest()
source("R/preprocessing.R")      # filter_near_zero_variance(), standardizer
source("R/model_smo_svm.R")      # train_smo_svm(), predict_smo_svm()
source("R/metrics.R")            # compute_binary_metrics()

set.seed(SEED)

# Pilot configuration -----------------------------------------------------
ACCESSION    <- "GSE25055"   # PILOT: discovery cohort only.
LABEL_FIELD  <- "pathologic_response_pcr_rd"
POSITIVE     <- "pCR"        # positive (minority) class
TOP_K        <- 100          # first pilot feature budget
N_FOLDS      <- 5            # cross-validation folds
N_REPEATS    <- 5            # repeated CV for a more stable estimate
COST         <- 1            # fixed linear-SVM cost (no tuning in leaky baseline)

RESULTS_DIR  <- file.path("results", "pilot_gse25055")
TABLES_DIR   <- file.path("tables", "pilot_gse25055")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)

# ---- Data loading -----------------------------------------------------------
# Loads GSE25055 expression matrix + pCR/RD labels via GEOquery.
# Returns list(x = samples-x-genes matrix, y = factor labels, ids = GSM ids).
load_gse25055 <- function(accession = ACCESSION, label_field = LABEL_FIELD) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package GEOquery (Bioconductor) is required to load ", accession, ".")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package Biobase (Bioconductor) is required.")
  }

  # getGEO returns a list of ExpressionSets (one per platform). The pilot uses
  # the single GPL96 series matrix; we take the first element.
  gset <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)
  eset <- if (length(gset) > 1) gset[[1]] else gset[[1]]

  expr <- Biobase::exprs(eset)        # genes (rows) x samples (cols)
  pheno <- Biobase::pData(eset)

  # --- Locate the pCR/RD label column robustly --------------------------------
  # GEO 'characteristics_ch1' fields are often split into columns whose values
  # look like "pathologic_response_pcr_rd: pCR". We search for the label field.
  label_col <- NULL
  for (cn in colnames(pheno)) {
    vals <- as.character(pheno[[cn]])
    if (any(grepl(label_field, vals, ignore.case = TRUE)) ||
        grepl(label_field, cn, ignore.case = TRUE)) {
      label_col <- cn
      break
    }
  }
  if (is.null(label_col)) {
    stop("Could not find label field '", label_field, "' in ", accession,
         " phenotype data. Inspect pData(eset) column names and values.")
  }

  raw_labels <- as.character(pheno[[label_col]])
  # Strip a leading "pathologic_response_pcr_rd: " prefix if present.
  raw_labels <- sub(paste0("^.*", label_field, "\\s*[:=]?\\s*"), "", raw_labels,
                    ignore.case = TRUE)
  raw_labels <- trimws(raw_labels)

  # --- Exclude NA-coded samples ----------------------------------------------
  # NA may appear as the string "NA", true NA, or empty. All are excluded.
  is_na <- is.na(raw_labels) | raw_labels %in% c("NA", "na", "N/A", "", "NaN")
  keep  <- !is_na

  message(sprintf("[load] %s: %d samples total; excluding %d NA-coded; keeping %d.",
                  accession, length(raw_labels), sum(is_na), sum(keep)))

  labels <- factor(raw_labels[keep], levels = c("RD", "pCR"))
  if (any(is.na(labels))) {
    stop("Unexpected label values after NA exclusion: ",
         paste(unique(raw_labels[keep]), collapse = ", "),
         ". Expected only 'pCR' / 'RD'. Inspect coding before proceeding.")
  }

  # Samples x genes matrix (rows = samples). Drop NA samples.
  x <- t(expr[, keep, drop = FALSE])
  ids <- colnames(expr)[keep]
  rownames(x) <- ids
  ord <- order(ids)
  x <- x[ord, , drop = FALSE]
  labels <- labels[ord]
  ids <- ids[ord]

  list(x = x, y = labels, ids = ids, label_col = label_col)
}

# ---- Light, label-independent preprocessing --------------------------------
# Removing constant/near-zero-variance probes does not use the outcome, so it
# is not a source of label leakage. (Feature SELECTION below IS leaky.)
prep_expression <- function(x) {
  nzv <- filter_near_zero_variance(x, cutoff = 1e-8)
  nzv$x
}

# ---- LEAKY step: global supervised feature selection ------------------------
# !!! LEAKAGE ON PURPOSE !!!
# Selecting features from the FULL dataset (all samples, all labels) means the
# CV folds below are evaluated on probes that were chosen with knowledge of the
# test-fold labels. This is the exact mistake the manuscript audits.
leaky_global_feature_selection <- function(x, y, top_k = TOP_K) {
  message("[LEAKY] Selecting top ", top_k,
          " features on the FULL dataset BEFORE cross-validation.")
  fs <- select_features(x, y, method = "t_test", top_k = top_k)
  fs$features
}

# ---- Cross-validation over the (leaked) fixed feature set -------------------
run_leaky_cv <- function(x, y, selected_features,
                         n_folds = N_FOLDS, n_repeats = N_REPEATS,
                         cost = COST, positive = POSITIVE) {
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package caret is required for fold creation.")
  }
  x_sel <- x[, selected_features, drop = FALSE]

  preds <- list()
  for (rep in seq_len(n_repeats)) {
    set.seed(SEED + rep)
    folds <- caret::createFolds(y, k = n_folds, list = TRUE, returnTrain = FALSE)
    for (f in seq_along(folds)) {
      test_idx  <- folds[[f]]
      train_idx <- setdiff(seq_along(y), test_idx)

      # NOTE: scaling here is still fit on the training rows only, but because
      # the FEATURE SET was chosen globally, the pipeline as a whole is leaky.
      std <- fit_standardizer(x_sel[train_idx, , drop = FALSE])
      x_tr <- apply_standardizer(x_sel[train_idx, , drop = FALSE], std)
      x_te <- apply_standardizer(x_sel[test_idx,  , drop = FALSE], std)

      model <- train_smo_svm(x_tr, y[train_idx], kernel = "linear",
                             cost = cost, class_weights = TRUE,
                             probability = TRUE)
      pr <- predict_smo_svm(model, x_te, positive_class = positive)

      preds[[length(preds) + 1]] <- data.frame(
        repeat_id = rep,
        fold      = f,
        sample_id = rownames(x_sel)[test_idx],
        truth     = as.character(y[test_idx]),
        predicted = pr$predicted_class,
        prob_pos  = pr$probability_positive,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, preds)
}

# ---- Summarize predictions into metrics -------------------------------------
summarize_metrics <- function(pred_df, positive = POSITIVE) {
  per_rep <- lapply(split(pred_df, pred_df$repeat_id), function(d) {
    m <- compute_binary_metrics(truth = d$truth, estimate = d$predicted,
                                probability = d$prob_pos, positive = positive)
    sens <- sum(d$truth == positive & d$predicted == positive) /
            sum(d$truth == positive)
    spec <- sum(d$truth != positive & d$predicted != positive) /
            sum(d$truth != positive)
    m$sensitivity <- sens
    m$specificity <- spec
    m
  })
  agg <- do.call(rbind, per_rep)
  data.frame(
    pipeline          = "A_leaky_baseline",
    auroc             = mean(agg$auroc, na.rm = TRUE),
    pr_auc            = mean(agg$pr_auc, na.rm = TRUE),
    balanced_accuracy = mean(agg$balanced_accuracy, na.rm = TRUE),
    mcc               = mean(agg$mcc, na.rm = TRUE),
    sensitivity       = mean(agg$sensitivity, na.rm = TRUE),
    specificity       = mean(agg$specificity, na.rm = TRUE),
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

  selected <- leaky_global_feature_selection(x, y, top_k = TOP_K)  # LEAKY
  pred_df  <- run_leaky_cv(x, y, selected)
  metrics  <- summarize_metrics(pred_df)

  # Persist pilot outputs (Pipeline A portion).
  utils::write.csv(metrics,
                   file.path(RESULTS_DIR, "leaky_baseline_metrics.csv"),
                   row.names = FALSE)
  utils::write.csv(pred_df,
                   file.path(RESULTS_DIR, "leaky_baseline_predictions.csv"),
                   row.names = FALSE)
  utils::write.csv(data.frame(rank = seq_along(selected), feature = selected),
                   file.path(RESULTS_DIR, "leaky_baseline_selected_features.csv"),
                   row.names = FALSE)

  message("[main] Leaky baseline complete. Metrics:")
  print(metrics)
  invisible(list(metrics = metrics, predictions = pred_df, features = selected))
}

if (sys.nframe() == 0) {
  main()
}
