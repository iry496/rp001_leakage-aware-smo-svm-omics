#!/usr/bin/env Rscript
# =============================================================================
# scripts/09_repeated_nested_cv.R
# -----------------------------------------------------------------------------
# Phase 2B: REPEATED NESTED CROSS-VALIDATION robustness check (GSE25055 ONLY).
#
# Tests whether the leakage gap (leaky - guarded) and the feature-stability
# results are robust across different random seeds / fold assignments, rather
# than artefacts of a single split. For each seed we run BOTH pipelines on the
# TRUE labels and record per-seed metrics, the leakage gap, and guarded-arm
# feature-stability metrics.
#
# This is a MATCHED SEED-LEVEL COMPARISON, NOT an identical fold structure:
# the leaky arm uses 5x5 repeated CV and the guarded arm uses nested 5-outer x
# 5-inner CV. They do not share an identical fold layout; instead, within each
# seed BOTH arms draw their stratified partitions from the same master seed, so
# the per-seed gap is a matched (paired) comparison under a common seed.
#
# Design (approved):
#   * GSE25055 only; true labels (no permutation).
#   * Leaky arm: global t-test top-K=100 BEFORE CV, then 5x5 repeated CV.
#   * Guarded arm: nested CV (5 outer x 5 inner); FS + tuning inside train folds.
#   * Feature selection: t-test top-K = 100 (unchanged). Cost grid 0.25/1/4.
#   * Vectorized Welch selector reused from script 08 (validated bit-identical),
#     with original-selector fallback toggle + a selector-validation guard.
#   * Parallelism ACROSS SEEDS ONLY (mclapply); CV loops stay serial.
#   * Anchor seed 20260620 must reproduce the committed pilot point estimates
#     and committed feature-stability numbers (hard gate).
#
# Scope guard: do NOT load GSE25065 / GSE41998 / GSE20194 / GSE20271.
# No manuscript / external / raw / processed data is committed by this script.
#
# Run (auto-runs main() via the sys.nframe guard; tag/n_seeds from args):
#   Rscript scripts/09_repeated_nested_cv.R smoke 3      # 3-seed smoke (default)
#   Rscript scripts/09_repeated_nested_cv.R full 30      # later, after review
# =============================================================================

suppressWarnings(suppressMessages({
  source("R/00_config.R")          # SEED (20260620), %||%
  source("R/feature_selection.R")  # select_features() (original t.test selector)
  source("R/preprocessing.R")      # filter_near_zero_variance(), standardizer
  source("R/model_smo_svm.R")      # train_smo_svm(), predict_smo_svm()
  source("R/metrics.R")            # compute_binary_metrics()  (pROC + PRROC)
  source("R/leakage_checks.R")     # assert_no_overlap()
}))

# ---- Run configuration ------------------------------------------------------
N_SEEDS           <- 3             # default = smoke (3 seeds incl. anchor)
RUN_TAG           <- "smoke"       # output filename tag (notes/figures)
USE_FAST_SELECTOR <- TRUE          # FALSE -> original t.test() selector
USE_PARALLEL      <- TRUE          # FALSE -> serial fallback (debug)
N_WORKERS         <- max(1, min(4, parallel::detectCores() - 2))

# Optional CLI overrides:
#   Rscript 09_repeated_nested_cv.R <RUN_TAG> <N_SEEDS>   # full run
#   Rscript 09_repeated_nested_cv.R posthoc <RUN_TAG>     # tests-only, NO CV rerun
.args <- commandArgs(trailingOnly = TRUE)
MODE <- "run"
if (length(.args) >= 1 && identical(.args[1], "posthoc")) {
  MODE <- "posthoc"
  if (length(.args) >= 2 && nzchar(.args[2])) RUN_TAG <- .args[2] else RUN_TAG <- "full"
} else {
  if (length(.args) >= 1 && nzchar(.args[1])) RUN_TAG <- .args[1]
  if (length(.args) >= 2 && nzchar(.args[2])) N_SEEDS <- as.integer(.args[2])
}

# Seed list: anchor SEED first, then contiguous. N_SEEDS=3 -> 20260620/21/22.
SEEDS <- SEED + seq_len(N_SEEDS) - 1L

# ---- Pilot configuration (identical to scripts 02 and 03) -------------------
ACCESSION   <- "GSE25055"
LABEL_FIELD <- "pathologic_response_pcr_rd"
POSITIVE    <- "pCR"
TOP_K       <- 100
KERNEL      <- "linear"
N_FOLDS     <- 5; N_REPEATS <- 5; COST <- 1                     # leaky arm
OUTER_FOLDS <- 5; INNER_FOLDS <- 5; COST_GRID <- c(0.25, 1, 4)  # guarded arm

# Committed references for the anchor reproduction gate (from Phase 1/2A):
REF_LEAKY_AUROC   <- 0.7830
REF_GUARDED_AUROC <- 0.7032
REF_NOGUEIRA      <- 0.5128
REF_JACCARD       <- 0.3487
REF_CORE          <- 26L
REF_TAIL          <- 105L

RESULTS_DIR <- file.path("results", "repeated_cv")
TABLES_DIR  <- file.path("tables",  "repeated_cv")
FIG_DIR     <- file.path("figures", "repeated_cv")
CACHE_X <- file.path("processed_data", "gse25055_perm_cache.rds")   # gitignored (shared w/ 08)
CACHE_Y <- file.path("processed_data", "gse25055_perm_labels.rds")  # gitignored
for (d in c(RESULTS_DIR, TABLES_DIR, FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Data loader (verbatim from scripts 03 / 08) ----------------------------
load_gse25055 <- function(accession = ACCESSION, label_field = LABEL_FIELD) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery required.")
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase required.")
  gset  <- GEOquery::getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)
  eset  <- gset[[1]]
  expr  <- Biobase::exprs(eset)
  pheno <- Biobase::pData(eset)
  label_col <- NULL
  for (cn in colnames(pheno)) {
    vals <- as.character(pheno[[cn]])
    if (any(grepl(label_field, vals, ignore.case = TRUE)) ||
        grepl(label_field, cn, ignore.case = TRUE)) { label_col <- cn; break }
  }
  if (is.null(label_col)) stop("Label field '", label_field, "' not found.")
  raw_labels <- as.character(pheno[[label_col]])
  raw_labels <- sub(paste0("^.*", label_field, "\\s*[:=]?\\s*"), "", raw_labels, ignore.case = TRUE)
  raw_labels <- trimws(raw_labels)
  is_na <- is.na(raw_labels) | raw_labels %in% c("NA", "na", "N/A", "", "NaN")
  keep  <- !is_na
  message(sprintf("[load] %s: %d total; excluding %d NA; keeping %d.",
                  accession, length(raw_labels), sum(is_na), sum(keep)))
  labels <- factor(raw_labels[keep], levels = c("RD", "pCR"))
  if (any(is.na(labels))) stop("Unexpected label values after NA exclusion.")
  x <- t(expr[, keep, drop = FALSE]); rownames(x) <- colnames(expr)[keep]
  ord <- order(rownames(x)); x <- x[ord, , drop = FALSE]; labels <- labels[ord]
  list(x = x, y = labels)
}

prep_expression <- function(x) filter_near_zero_variance(x, cutoff = 1e-8)$x

load_or_cache <- function() {
  if (file.exists(CACHE_X) && file.exists(CACHE_Y)) {
    message("[cache] reading processed matrix/labels from rds cache.")
    x <- readRDS(CACHE_X); y <- readRDS(CACHE_Y)
    ord <- order(rownames(x))
    return(list(x = x[ord, , drop = FALSE], y = y[ord]))
  }
  dat <- load_gse25055(); x <- prep_expression(dat$x); y <- dat$y
  dir.create("processed_data", recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, CACHE_X); saveRDS(y, CACHE_Y)
  message("[cache] wrote rds cache (gitignored).")
  list(x = x, y = y)
}

# ---- Feature selectors (verbatim from script 08) ----------------------------
select_ttest_fast <- function(x_train, y_train, top_k = TOP_K) {
  x <- as.matrix(x_train); y <- as.factor(y_train)
  lv <- levels(y); g1 <- which(y == lv[1]); g2 <- which(y == lv[2])
  n1 <- length(g1); n2 <- length(g2)
  m1 <- colMeans(x[g1, , drop = FALSE]); m2 <- colMeans(x[g2, , drop = FALSE])
  v1 <- (colSums(x[g1, , drop = FALSE]^2) - n1 * m1^2) / (n1 - 1)
  v2 <- (colSums(x[g2, , drop = FALSE]^2) - n2 * m2^2) / (n2 - 1)
  v1[v1 < 0] <- 0; v2[v2 < 0] <- 0
  se <- sqrt(v1 / n1 + v2 / n2)
  scores <- abs((m1 - m2) / se)
  scores[!is.finite(scores)] <- 0
  names(scores) <- colnames(x)
  selected <- names(sort(scores, decreasing = TRUE))[seq_len(min(top_k, length(scores)))]
  list(features = selected, scores = scores[selected], method = "t_test_fast")
}

fs_select <- function(x_train, y_train, top_k = TOP_K) {
  if (isTRUE(USE_FAST_SELECTOR)) select_ttest_fast(x_train, y_train, top_k)$features
  else select_features(x_train, y_train, method = "t_test", top_k = top_k)$features
}

validate_selector <- function(x, y, n_check = 5) {
  cases <- c(list(identity = y),
             setNames(lapply(seq_len(n_check), function(b) {
               set.seed(SEED + 900000 + b); y[sample.int(length(y))]
             }), paste0("perm", seq_len(n_check))))
  do.call(rbind, lapply(names(cases), function(nm) {
    yy  <- cases[[nm]]
    old <- select_features(x, yy, method = "t_test", top_k = TOP_K)$features
    new <- select_ttest_fast(x, yy, top_k = TOP_K)$features
    data.frame(case = nm, top_k = TOP_K,
               overlap = length(intersect(old, new)),
               exact_order = identical(old, new), stringsAsFactors = FALSE)
  }))
}

# ---- LEAKY arm (logic from script 02/08; seed-parameterized) ----------------
run_leaky_cv <- function(x, y, selected_features, seed) {
  x_sel <- x[, selected_features, drop = FALSE]; preds <- list()
  for (rep in seq_len(N_REPEATS)) {
    set.seed(seed + rep)
    folds <- caret::createFolds(y, k = N_FOLDS, list = TRUE, returnTrain = FALSE)
    for (f in seq_along(folds)) {
      test_idx <- folds[[f]]; train_idx <- setdiff(seq_along(y), test_idx)
      std  <- fit_standardizer(x_sel[train_idx, , drop = FALSE])
      x_tr <- apply_standardizer(x_sel[train_idx, , drop = FALSE], std)
      x_te <- apply_standardizer(x_sel[test_idx,  , drop = FALSE], std)
      model <- train_smo_svm(x_tr, y[train_idx], kernel = KERNEL, cost = COST,
                             class_weights = TRUE, probability = TRUE)
      pr <- predict_smo_svm(model, x_te, positive_class = POSITIVE)
      preds[[length(preds) + 1]] <- data.frame(
        repeat_id = rep, fold = f, sample_id = rownames(x_sel)[test_idx],
        truth = as.character(y[test_idx]), predicted = pr$predicted_class,
        prob_pos = pr$probability_positive, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, preds)
}

leaky_summarize <- function(pred_df) {
  per_rep <- lapply(split(pred_df, pred_df$repeat_id), function(d) {
    m <- compute_binary_metrics(d$truth, d$predicted, d$prob_pos, POSITIVE)
    m$sensitivity <- sum(d$truth == POSITIVE & d$predicted == POSITIVE) / sum(d$truth == POSITIVE)
    m$specificity <- sum(d$truth != POSITIVE & d$predicted != POSITIVE) / sum(d$truth != POSITIVE)
    m
  })
  agg <- do.call(rbind, per_rep)
  c(auroc = mean(agg$auroc), pr_auc = mean(agg$pr_auc),
    balanced_accuracy = mean(agg$balanced_accuracy), mcc = mean(agg$mcc),
    sensitivity = mean(agg$sensitivity), specificity = mean(agg$specificity))
}

# ---- GUARDED nested arm (logic from script 03/08; seed-parameterized) --------
tune_inner <- function(x_train, y_train, seed) {
  set.seed(seed)
  folds <- caret::createFolds(y_train, k = INNER_FOLDS, returnTrain = FALSE)
  score_for_cost <- function(cost) {
    bal <- numeric(0)
    for (f in seq_along(folds)) {
      val_idx <- folds[[f]]; tr_idx <- setdiff(seq_along(y_train), val_idx)
      feats <- fs_select(x_train[tr_idx, , drop = FALSE], y_train[tr_idx], TOP_K)
      std  <- fit_standardizer(x_train[tr_idx, feats, drop = FALSE])
      x_tr <- apply_standardizer(x_train[tr_idx,  feats, drop = FALSE], std)
      x_va <- apply_standardizer(x_train[val_idx, feats, drop = FALSE], std)
      model <- train_smo_svm(x_tr, y_train[tr_idx], kernel = KERNEL, cost = cost,
                             class_weights = TRUE, probability = TRUE)
      pr <- predict_smo_svm(model, x_va, positive_class = POSITIVE)
      sens <- sum(y_train[val_idx] == POSITIVE & pr$predicted_class == POSITIVE) / sum(y_train[val_idx] == POSITIVE)
      spec <- sum(y_train[val_idx] != POSITIVE & pr$predicted_class != POSITIVE) / sum(y_train[val_idx] != POSITIVE)
      bal <- c(bal, mean(c(sens, spec), na.rm = TRUE))
    }
    mean(bal, na.rm = TRUE)
  }
  COST_GRID[which.max(vapply(COST_GRID, score_for_cost, numeric(1)))]
}

run_nested <- function(x, y, seed) {
  set.seed(seed)
  folds <- caret::createFolds(y, k = OUTER_FOLDS, returnTrain = FALSE)
  preds <- list(); fold_feats <- vector("list", length(folds))
  assign_rows <- list(); best_costs <- numeric(length(folds))
  for (f in seq_along(folds)) {
    test_idx <- folds[[f]]; train_idx <- setdiff(seq_along(y), test_idx)
    assert_no_overlap(rownames(x)[train_idx], rownames(x)[test_idx])
    x_train <- x[train_idx, , drop = FALSE]; y_train <- y[train_idx]
    best_cost <- tune_inner(x_train, y_train, seed)
    feats <- fs_select(x_train, y_train, TOP_K)
    fold_feats[[f]] <- feats; best_costs[f] <- best_cost
    std   <- fit_standardizer(x_train[, feats, drop = FALSE])
    x_tr  <- apply_standardizer(x_train[, feats, drop = FALSE], std)
    x_te  <- apply_standardizer(x[test_idx, feats, drop = FALSE], std)
    model <- train_smo_svm(x_tr, y_train, kernel = KERNEL, cost = best_cost,
                           class_weights = TRUE, probability = TRUE)
    pr <- predict_smo_svm(model, x_te, positive_class = POSITIVE)
    preds[[f]] <- data.frame(
      fold = f, sample_id = rownames(x)[test_idx], truth = as.character(y[test_idx]),
      predicted = pr$predicted_class, prob_pos = pr$probability_positive,
      stringsAsFactors = FALSE)
    assign_rows[[f]] <- data.frame(fold = f, sample_id = rownames(x)[test_idx],
                                   stringsAsFactors = FALSE)
  }
  list(predictions = do.call(rbind, preds), fold_feats = fold_feats,
       fold_assign = do.call(rbind, assign_rows), best_costs = best_costs)
}

guarded_summarize <- function(pred_df) {
  m <- compute_binary_metrics(pred_df$truth, pred_df$predicted, pred_df$prob_pos, POSITIVE)
  c(auroc = m$auroc, pr_auc = m$pr_auc, balanced_accuracy = m$balanced_accuracy, mcc = m$mcc,
    sensitivity = sum(pred_df$truth == POSITIVE & pred_df$predicted == POSITIVE) / sum(pred_df$truth == POSITIVE),
    specificity = sum(pred_df$truth != POSITIVE & pred_df$predicted != POSITIVE) / sum(pred_df$truth != POSITIVE))
}

# ---- feature-stability metrics (formulas from script 04) --------------------
stability_from_feats <- function(fold_feats, P) {
  M <- length(fold_feats)
  k_each <- vapply(fold_feats, length, integer(1)); kbar <- mean(k_each)
  freq_int <- as.integer(table(unlist(fold_feats)))
  phat      <- freq_int / M
  var_terms <- (M / (M - 1)) * phat * (1 - phat)
  nogueira  <- 1 - (sum(var_terms) / P) / ((kbar / P) * (1 - kbar / P))
  pair <- numeric(0)
  for (i in 1:(M - 1)) for (j in (i + 1):M) {
    a <- fold_feats[[i]]; b <- fold_feats[[j]]
    pair <- c(pair, length(intersect(a, b)) / length(union(a, b)))
  }
  c(nogueira = nogueira, mean_jaccard = mean(pair), median_jaccard = median(pair),
    total_unique = length(freq_int),
    stable_core_count = sum(freq_int == M), unstable_tail_count = sum(freq_int == 1))
}

# ---- one seed (both arms + stability) -> list of result pieces ---------------
seed_worker <- function(seed, x, y, P) {
  tstart <- Sys.time(); is_anchor <- (seed == SEED)
  # leaky arm (global FS is seed-independent; CV fold assignment varies by seed)
  lk_pred <- run_leaky_cv(x, y, fs_select(x, y, TOP_K), seed)
  lk <- leaky_summarize(lk_pred)
  # guarded nested arm
  nest <- run_nested(x, y, seed)
  gd   <- guarded_summarize(nest$predictions)
  st   <- stability_from_feats(nest$fold_feats, P)
  runtime <- as.numeric(difftime(Sys.time(), tstart, units = "secs"))

  gap_row <- data.frame(
    seed = seed, is_anchor = is_anchor,
    leaky_auroc = lk["auroc"], leaky_pr_auc = lk["pr_auc"],
    leaky_balanced_accuracy = lk["balanced_accuracy"], leaky_mcc = lk["mcc"],
    leaky_sensitivity = lk["sensitivity"], leaky_specificity = lk["specificity"],
    guarded_auroc = gd["auroc"], guarded_pr_auc = gd["pr_auc"],
    guarded_balanced_accuracy = gd["balanced_accuracy"], guarded_mcc = gd["mcc"],
    guarded_sensitivity = gd["sensitivity"], guarded_specificity = gd["specificity"],
    gap_auroc = lk["auroc"] - gd["auroc"], gap_pr_auc = lk["pr_auc"] - gd["pr_auc"],
    guarded_best_costs = paste(nest$best_costs, collapse = "|"),
    runtime_secs = runtime, stringsAsFactors = FALSE)

  stab_row <- data.frame(
    seed = seed, is_anchor = is_anchor,
    nogueira_stability = st["nogueira"], mean_jaccard = st["mean_jaccard"],
    median_jaccard = st["median_jaccard"], total_unique_features = st["total_unique"],
    stable_core_count = st["stable_core_count"], unstable_tail_count = st["unstable_tail_count"],
    top_k = TOP_K, n_outer_folds = OUTER_FOLDS, feature_universe = P,
    stringsAsFactors = FALSE)

  feat_rows <- do.call(rbind, lapply(seq_along(nest$fold_feats), function(f)
    data.frame(seed = seed, fold = f, rank = seq_along(nest$fold_feats[[f]]),
               feature = nest$fold_feats[[f]], stringsAsFactors = FALSE)))

  leaky_assign <- data.frame(seed = seed, arm = "leaky",
                             repeat_id = lk_pred$repeat_id, fold = lk_pred$fold,
                             sample_id = lk_pred$sample_id, stringsAsFactors = FALSE)
  guarded_assign <- data.frame(seed = seed, arm = "guarded_outer",
                               repeat_id = NA_integer_, fold = nest$fold_assign$fold,
                               sample_id = nest$fold_assign$sample_id, stringsAsFactors = FALSE)

  list(gap = gap_row, stab = stab_row, feats = feat_rows,
       folds = rbind(leaky_assign, guarded_assign))
}

# ---- figures (robust to small n) --------------------------------------------
write_figures <- function(gap, stab) {
  idx <- seq_len(nrow(gap))
  for (dev_fun in list(
    function() png(file.path(FIG_DIR, sprintf("leakage_gap_distribution_%s.png", RUN_TAG)), width = 1100, height = 500),
    function() pdf(file.path(FIG_DIR, sprintf("leakage_gap_distribution_%s.pdf", RUN_TAG)), width = 11, height = 5))) {
    dev_fun(); par(mfrow = c(1, 2))
    # Panel 1: per-seed leaky vs guarded AUROC
    yr <- range(c(gap$leaky_auroc, gap$guarded_auroc, REF_LEAKY_AUROC, REF_GUARDED_AUROC))
    plot(idx, gap$leaky_auroc, pch = 19, col = "#C0392B", ylim = yr, xaxt = "n",
         xlab = "seed index", ylab = "AUROC",
         main = sprintf("Per-seed AUROC (%d seeds)", nrow(gap)))
    points(idx, gap$guarded_auroc, pch = 17, col = "#1E8449")
    axis(1, at = idx, labels = idx)
    abline(h = REF_LEAKY_AUROC, lty = 3, col = "#C0392B")
    abline(h = REF_GUARDED_AUROC, lty = 3, col = "#1E8449")
    legend("right", bty = "n", pch = c(19, 17), col = c("#C0392B", "#1E8449"),
           legend = c("leaky", "guarded"))
    # Panel 2: per-seed leakage gap (boxplot when n>=10, else points)
    if (nrow(gap) >= 10) {
      boxplot(gap$gap_auroc, col = "#D9E2F3", ylab = "delta AUROC (leaky - guarded)",
              main = "Leakage gap across seeds")
      points(rep(1, nrow(gap)), gap$gap_auroc, pch = 19, col = "#34495E")
    } else {
      plot(idx, gap$gap_auroc, pch = 19, col = "#34495E", xaxt = "n",
           xlab = "seed index", ylab = "delta AUROC (leaky - guarded)",
           main = "Leakage gap across seeds")
      axis(1, at = idx, labels = idx)
    }
    abline(h = 0, lty = 3)
    dev.off()
  }
  for (dev_fun in list(
    function() png(file.path(FIG_DIR, sprintf("stability_distribution_%s.png", RUN_TAG)), width = 1100, height = 500),
    function() pdf(file.path(FIG_DIR, sprintf("stability_distribution_%s.pdf", RUN_TAG)), width = 11, height = 5))) {
    dev_fun(); par(mfrow = c(1, 2))
    plot(idx, stab$nogueira_stability, pch = 19, col = "#6C3483", xaxt = "n",
         ylim = range(c(stab$nogueira_stability, REF_NOGUEIRA)),
         xlab = "seed index", ylab = "Nogueira stability",
         main = "Per-seed Nogueira stability")
    axis(1, at = idx, labels = idx); abline(h = REF_NOGUEIRA, lty = 3, col = "#6C3483")
    plot(idx, stab$mean_jaccard, pch = 19, col = "#1F618D", xaxt = "n",
         ylim = range(c(stab$mean_jaccard, REF_JACCARD)),
         xlab = "seed index", ylab = "mean pairwise Jaccard",
         main = "Per-seed mean Jaccard")
    axis(1, at = idx, labels = idx); abline(h = REF_JACCARD, lty = 3, col = "#1F618D")
    dev.off()
  }
}

# ---- gap significance: Wilcoxon signed-rank (paired delta vs 0) -------------
# Canonical test for the leakage gap. One-sample Wilcoxon on the per-seed
# deltas against mu = 0. exact = FALSE (ties present at n = 30) -> normal
# approximation with continuity correction (base stats::wilcox.test).
gap_tests_df <- function(gap) {
  mk <- function(label, v) {
    tt <- suppressWarnings(stats::wilcox.test(v, mu = 0, alternative = "two.sided",
                                              exact = FALSE, correct = TRUE))
    data.frame(
      statistic   = label,
      n_seeds     = length(v),
      n_nonzero   = sum(v != 0),
      n_positive  = sum(v > 0),
      n_negative  = sum(v < 0),
      median      = round(median(v), 4),
      ci2.5       = round(unname(quantile(v, 0.025)), 4),
      ci97.5      = round(unname(quantile(v, 0.975)), 4),
      V_statistic = unname(tt$statistic),
      p_value     = signif(tt$p.value, 6),
      alternative = "two.sided",
      method      = "Wilcoxon signed-rank (normal approx, continuity corrected)",
      stringsAsFactors = FALSE)
  }
  rbind(mk("delta_auroc", gap$gap_auroc), mk("delta_pr_auc", gap$gap_pr_auc))
}

gap_tests_notes_lines <- function(gap, tests) {
  ta <- tests[tests$statistic == "delta_auroc", ]
  tp <- tests[tests$statistic == "delta_pr_auc", ]
  c("", "## Gap significance (Wilcoxon signed-rank, paired delta vs 0)",
    sprintf("- delta AUROC: V = %.0f, p = %.3g; positive in %d/%d seeds.",
            ta$V_statistic, ta$p_value, ta$n_positive, ta$n_seeds),
    sprintf("- delta PR-AUC: V = %.0f, p = %.3g; positive in %d/%d seeds.",
            tp$V_statistic, tp$p_value, tp$n_positive, tp$n_seeds),
    sprintf(paste0("- The workflow contrasts are positive in %d/%d seeds for AUROC and %d/%d seeds for PR-AUC ",
                   "(Wilcoxon p = %.3g and %.3g, respectively), but the ",
                   "2.5-97.5%% interval includes small negative values, so claims should remain ",
                   "cautious: this is a within-cohort, methodology/diagnostic finding, not ",
                   "biomarker discovery."),
            ta$n_positive, ta$n_seeds, tp$n_positive, tp$n_seeds,
            ta$p_value, tp$p_value))
}

# ---- driver -----------------------------------------------------------------
main <- function() {
  t0 <- Sys.time()
  message(sprintf("[cfg] RUN_TAG=%s N_SEEDS=%d FAST=%s PARALLEL=%s N_WORKERS=%d",
                  RUN_TAG, N_SEEDS, USE_FAST_SELECTOR, USE_PARALLEL, N_WORKERS))
  message(sprintf("[cfg] seeds: %s", paste(SEEDS, collapse = ", ")))
  dd <- load_or_cache(); x <- dd$x; y <- dd$y; P <- ncol(x)
  message(sprintf("[rcv] N=%d; P=%d; class balance %s", length(y), P,
                  paste(names(table(y)), table(y), sep = "=", collapse = ", ")))

  # selector validation (only meaningful when the fast selector is used)
  val <- NULL
  if (isTRUE(USE_FAST_SELECTOR)) {
    message("[validate] comparing original vs fast selector (identity + 5 shuffles)...")
    val <- validate_selector(x, y, n_check = 5); print(val)
  }

  run_one <- function(s) seed_worker(s, x, y, P)
  if (isTRUE(USE_PARALLEL) && .Platform$OS.type != "windows" && N_WORKERS > 1) {
    message(sprintf("[run] parallel across seeds (mclapply, %d workers); CV loops serial.", N_WORKERS))
    pieces <- parallel::mclapply(SEEDS, run_one, mc.cores = N_WORKERS, mc.preschedule = FALSE)
  } else {
    message("[run] serial across seeds."); pieces <- lapply(SEEDS, run_one)
  }

  gap  <- do.call(rbind, lapply(pieces, `[[`, "gap"))
  stab <- do.call(rbind, lapply(pieces, `[[`, "stab"))
  feats <- do.call(rbind, lapply(pieces, `[[`, "feats"))
  folds <- do.call(rbind, lapply(pieces, `[[`, "folds"))
  gap  <- gap[order(gap$seed), ];  rownames(gap)  <- NULL
  stab <- stab[order(stab$seed), ]; rownames(stab) <- NULL

  write.csv(gap,  file.path(TABLES_DIR,  "leakage_gap_by_seed.csv"), row.names = FALSE)
  write.csv(stab, file.path(TABLES_DIR,  "stability_by_seed.csv"),   row.names = FALSE)
  write.csv(feats, file.path(RESULTS_DIR, "selected_features_by_seed_fold.csv"), row.names = FALSE)
  write.csv(folds, file.path(RESULTS_DIR, "fold_assignments_by_seed.csv"),       row.names = FALSE)
  write_figures(gap, stab)

  # canonical gap significance tests (Wilcoxon signed-rank, paired delta vs 0)
  tests <- gap_tests_df(gap)
  write.csv(tests, file.path(RESULTS_DIR, "repeated_cv_gap_tests.csv"), row.names = FALSE)

  # ---- anchor reproduction gate (seed == SEED) ------------------------------
  a_gap  <- gap[gap$seed == SEED, ]; a_stab <- stab[stab$seed == SEED, ]
  chk <- c(
    leaky_auroc   = abs(a_gap$leaky_auroc   - REF_LEAKY_AUROC)   < 0.01,
    guarded_auroc = abs(a_gap$guarded_auroc - REF_GUARDED_AUROC) < 0.01,
    nogueira      = abs(a_stab$nogueira_stability - REF_NOGUEIRA) < 0.01,
    mean_jaccard  = abs(a_stab$mean_jaccard - REF_JACCARD) < 0.01,
    stable_core   = a_stab$stable_core_count   == REF_CORE,
    unstable_tail = a_stab$unstable_tail_count == REF_TAIL)
  repro_ok <- all(chk)

  # ---- notes ----------------------------------------------------------------
  mean_per_seed <- mean(gap$runtime_secs)
  total_wall    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  wall_per_seed <- total_wall / nrow(gap)
  est <- function(n) sprintf("%.1f min (%.2f h)", n * wall_per_seed / 60, n * wall_per_seed / 3600)
  fmt <- function(v) paste(sprintf("%.4f", v), collapse = ", ")

  notes <- c(
    sprintf("# Repeated Nested CV - %s run (GSE25055 only)", RUN_TAG), "",
    sprintf("- Seeds (%d): %s", nrow(gap), paste(gap$seed, collapse = ", ")),
    "- Matched seed-level comparison (NOT identical fold structure): leaky = 5x5 repeated CV,",
    "  guarded = nested 5-outer x 5-inner; within each seed both arms draw stratified folds",
    "  from the same master seed.",
    sprintf("- Fast Welch selector: %s ; parallel across seeds: %s (%d workers); CV loops serial.",
            USE_FAST_SELECTOR, USE_PARALLEL, N_WORKERS),
    "- Leaky arm: global t-test top-100 before CV. Guarded arm: t-test inside training folds only.",
    "- Cost grid: 0.25, 1, 4. Top-K = 100.",
    "", "## Anchor reproduction check (seed 20260620)",
    sprintf("- leaky AUROC = %.4f (ref %.4f)", a_gap$leaky_auroc, REF_LEAKY_AUROC),
    sprintf("- guarded AUROC = %.4f (ref %.4f)", a_gap$guarded_auroc, REF_GUARDED_AUROC),
    sprintf("- Nogueira = %.4f (ref %.4f)", a_stab$nogueira_stability, REF_NOGUEIRA),
    sprintf("- mean Jaccard = %.4f (ref %.4f)", a_stab$mean_jaccard, REF_JACCARD),
    sprintf("- stable-core count = %d (ref %d)", a_stab$stable_core_count, REF_CORE),
    sprintf("- unstable-tail count = %d (ref %d)", a_stab$unstable_tail_count, REF_TAIL),
    sprintf("- Overall: %s", if (repro_ok) "PASS" else "REVIEW (anchor did not reproduce within tolerance)"))
  if (!is.null(val)) notes <- c(notes, "", "## Selector validation (original vs fast Welch)",
    sprintf("- %s: overlap %d/%d, exact_order=%s.", val$case, val$overlap, val$top_k, val$exact_order))
  notes <- c(notes, "", "## Per-seed leakage gap (AUROC)",
    sprintf("- leaky AUROC: %s", fmt(gap$leaky_auroc)),
    sprintf("- guarded AUROC: %s", fmt(gap$guarded_auroc)),
    sprintf("- gap AUROC (leaky-guarded): %s", fmt(gap$gap_auroc)),
    sprintf("- gap PR-AUC: %s", fmt(gap$gap_pr_auc)),
    "", "## Per-seed feature stability",
    sprintf("- Nogueira: %s", fmt(stab$nogueira_stability)),
    sprintf("- mean Jaccard: %s", fmt(stab$mean_jaccard)),
    sprintf("- stable-core count: %s", paste(stab$stable_core_count, collapse = ", ")),
    sprintf("- unstable-tail count: %s", paste(stab$unstable_tail_count, collapse = ", ")),
    "", "## Runtime",
    sprintf("- mean %.1f s/seed (compute, within worker); wall %.1f s/seed; total %.1f min.",
            mean_per_seed, wall_per_seed, total_wall / 60),
    sprintf("- estimated full 20 seeds: %s ; 30 seeds: %s (at current wall throughput).", est(20), est(30)))
  notes <- c(notes, gap_tests_notes_lines(gap, tests))
  is_smoke <- grepl("smoke", RUN_TAG, ignore.case = TRUE) || nrow(gap) < 10
  footer <- if (is_smoke) {
    "_Smoke run: small number of seeds, intended for code validation and anchor reproduction only. Not the final robustness estimate._"
  } else {
    sprintf("_Full repeated-CV run across %d seeds; estimates seed/fold robustness of the leakage gap and feature stability (within-cohort, diagnostic; not biomarker discovery)._", nrow(gap))
  }
  notes <- c(notes, "", footer)
  writeLines(notes, file.path(RESULTS_DIR, sprintf("repeated_cv_%s_notes.md", RUN_TAG)))

  message(sprintf("[rcv] DONE (%s). Anchor reproduction: %s. mean %.1fs/seed; wall %.1fs/seed (%d workers). 20 ~ %s; 30 ~ %s.",
                  RUN_TAG, if (repro_ok) "PASS" else "REVIEW", mean_per_seed, wall_per_seed, N_WORKERS, est(20), est(30)))
  if (!repro_ok) {
    message("[rcv] WARNING: anchor reproduction FAILED. Failing checks: ",
            paste(names(chk)[!chk], collapse = ", "), ". Debug before the full run.")
  }
  invisible(list(gap = gap, stability = stab, selector_validation = val, repro_ok = repro_ok))
}

# ---- posthoc: compute gap tests + update notes WITHOUT rerunning the CV -----
# Reads the already-written leakage_gap_by_seed.csv, writes the canonical
# repeated_cv_gap_tests.csv, and adds the Wilcoxon section to the existing notes
# (also correcting a mislabeled smoke footer on a full run). Does NOT touch the
# leakage/stability CSV values.
main_posthoc <- function() {
  gap_path <- file.path(TABLES_DIR, "leakage_gap_by_seed.csv")
  if (!file.exists(gap_path)) stop("posthoc: ", gap_path, " not found; run the full CV first.")
  gap <- utils::read.csv(gap_path, stringsAsFactors = FALSE)
  tests <- gap_tests_df(gap)
  write.csv(tests, file.path(RESULTS_DIR, "repeated_cv_gap_tests.csv"), row.names = FALSE)

  notes_path <- file.path(RESULTS_DIR, sprintf("repeated_cv_%s_notes.md", RUN_TAG))
  if (file.exists(notes_path)) {
    lines <- readLines(notes_path)
    lines <- lines[!grepl("^_Smoke run:", lines)]                 # drop mislabeled footer
    lines <- lines[!grepl("^_Full repeated-CV run across", lines)] # drop any prior full footer
    while (length(lines) && !nzchar(lines[length(lines)])) lines <- lines[-length(lines)]
    if (!any(grepl("Gap significance", lines))) lines <- c(lines, gap_tests_notes_lines(gap, tests))
    lines <- c(lines, "", sprintf(
      "_Full repeated-CV run across %d seeds; estimates seed/fold robustness of the leakage gap and feature stability (within-cohort, diagnostic; not biomarker discovery)._",
      nrow(gap)))
    writeLines(lines, notes_path)
    message(sprintf("[posthoc] updated %s with Wilcoxon section.", basename(notes_path)))
  } else {
    message("[posthoc] notes file not found; wrote gap tests CSV only: ", notes_path)
  }
  message("[posthoc] repeated_cv_gap_tests.csv written (no CV rerun).")
  print(tests)
  invisible(tests)
}

if (sys.nframe() == 0) {
  if (identical(MODE, "posthoc")) main_posthoc() else main()
}
