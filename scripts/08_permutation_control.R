#!/usr/bin/env Rscript
# =============================================================================
# scripts/08_permutation_control.R
# -----------------------------------------------------------------------------
# Phase 2A: label-shuffle NEGATIVE CONTROL for the leakage audit (GSE25055 ONLY).
#
# Permutes the pCR/RD labels (sample-level) and reruns BOTH pipelines so that:
#   * the GUARDED nested pipeline should collapse to ~chance, and
#   * the LEAKY baseline should remain ABOVE chance even with random labels,
# supporting the interpretation that feature-selection leakage can inflate
# apparent performance (diagnostic only; not evidence of biological signal).
#
# Scientific design (UNCHANGED from the verified smoke run):
#   * GSE25055 only; label permutation only (expression matrix never altered).
#   * Leaky arm: global t-test top-K=100 on the FULL data BEFORE CV.
#   * Guarded arm: t-test top-K=100 inside training folds only.
#   * Same SVM model, same fold design (stratified, regenerated per permutation),
#     same metrics. Permutation 0 = identity (must reproduce naive ~0.7830,
#     guarded ~0.7032).
#
# Phase-2A optimizations (do not change the science):
#   1. Vectorized Welch t-test selector (USE_FAST_SELECTOR) — same statistic as
#      the pilot's t.test()-based selector, just computed by matrix algebra.
#   2. Selector validation: old vs fast on identity + 5 shuffles (overlap/order).
#   3. Local rds cache of the prepped matrix/labels (NOT committed; gitignored).
#   4. Parallelism ACROSS permutations only (USE_PARALLEL, N_WORKERS); CV loops
#      stay serial. RNG is reset inside the run functions, so results are
#      identical to serial.
#   5. Serial + original-selector fallbacks via toggles.
#
# Scope guard: do NOT load GSE25065 / GSE41998 / GSE20194 / GSE20271.
# No existing result / prediction / manuscript file is modified.
#
# Run (auto-runs main() via the sys.nframe guard; tag/B from args):
#   Rscript scripts/08_permutation_control.R smoke_fast_debug 5
#   Rscript scripts/08_permutation_control.R smoke_fast 20
#   Rscript scripts/08_permutation_control.R b200 200        # later, after review
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
B_PERM            <- 20            # shuffled permutations (default)
RUN_TAG           <- "smoke_fast"  # output filename tag (default)
USE_FAST_SELECTOR <- TRUE          # FALSE -> original t.test() selector
USE_PARALLEL      <- TRUE          # FALSE -> serial fallback (debug)
N_WORKERS         <- max(1, min(4, parallel::detectCores() - 2))

# Optional CLI overrides:  Rscript 08_permutation_control.R <RUN_TAG> <B_PERM>
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1 && nzchar(.args[1])) RUN_TAG <- .args[1]
if (length(.args) >= 2 && nzchar(.args[2])) B_PERM  <- as.integer(.args[2])

# ---- Pilot configuration (identical to scripts 02 and 03) -------------------
ACCESSION   <- "GSE25055"
LABEL_FIELD <- "pathologic_response_pcr_rd"
POSITIVE    <- "pCR"
TOP_K       <- 100
KERNEL      <- "linear"
N_FOLDS     <- 5; N_REPEATS <- 5; COST <- 1                 # leaky arm
OUTER_FOLDS <- 5; INNER_FOLDS <- 5; COST_GRID <- c(0.25, 1, 4)  # guarded arm
REF_LEAKY_AUROC   <- 0.7830
REF_GUARDED_AUROC <- 0.7032
PERM_SEED_OFFSET  <- 100000

RESULTS_DIR <- file.path("results", "permutation")
FIG_DIR     <- file.path("figures")
CACHE_X <- file.path("processed_data", "gse25055_perm_cache.rds")   # gitignored
CACHE_Y <- file.path("processed_data", "gse25055_perm_labels.rds")  # gitignored
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,     recursive = TRUE, showWarnings = FALSE)

# ---- Data loader (verbatim from script 03) ----------------------------------
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

# ---- Feature selectors ------------------------------------------------------
# Fast vectorized Welch two-sample t (same statistic as stats::t.test default).
select_ttest_fast <- function(x_train, y_train, top_k = TOP_K) {
  x <- as.matrix(x_train); y <- as.factor(y_train)
  lv <- levels(y); g1 <- which(y == lv[1]); g2 <- which(y == lv[2])
  n1 <- length(g1); n2 <- length(g2)
  m1 <- colMeans(x[g1, , drop = FALSE]); m2 <- colMeans(x[g2, , drop = FALSE])
  # per-group variance via sums of squares (vectorized, base R)
  v1 <- (colSums(x[g1, , drop = FALSE]^2) - n1 * m1^2) / (n1 - 1)
  v2 <- (colSums(x[g2, , drop = FALSE]^2) - n2 * m2^2) / (n2 - 1)
  v1[v1 < 0] <- 0; v2[v2 < 0] <- 0                  # guard fp negatives
  se <- sqrt(v1 / n1 + v2 / n2)
  tstat <- (m1 - m2) / se
  scores <- abs(tstat)
  scores[!is.finite(scores)] <- 0                   # NA/Inf/zero-var -> 0 (matches t.test catch)
  names(scores) <- colnames(x)
  selected <- names(sort(scores, decreasing = TRUE))[seq_len(min(top_k, length(scores)))]
  list(features = selected, scores = scores[selected], method = "t_test_fast")
}

# dispatcher used by both arms
fs_select <- function(x_train, y_train, top_k = TOP_K) {
  if (isTRUE(USE_FAST_SELECTOR)) select_ttest_fast(x_train, y_train, top_k)$features
  else select_features(x_train, y_train, method = "t_test", top_k = top_k)$features
}

# old-vs-fast selector validation (identity + n_check shuffles, global selection)
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

# ---- LEAKY arm (verbatim logic from script 02; FS routed via fs_select) ------
run_leaky_cv <- function(x, y, selected_features) {
  x_sel <- x[, selected_features, drop = FALSE]; preds <- list()
  for (rep in seq_len(N_REPEATS)) {
    set.seed(SEED + rep)
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
        repeat_id = rep, sample_id = rownames(x_sel)[test_idx],
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

run_leaky <- function(x, y) leaky_summarize(run_leaky_cv(x, y, fs_select(x, y, TOP_K)))

# ---- GUARDED nested arm (verbatim logic from script 03; FS via fs_select) ----
tune_inner <- function(x_train, y_train) {
  set.seed(SEED)
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

run_nested <- function(x, y) {
  set.seed(SEED)
  folds <- caret::createFolds(y, k = OUTER_FOLDS, returnTrain = FALSE)
  preds <- list()
  for (f in seq_along(folds)) {
    test_idx <- folds[[f]]; train_idx <- setdiff(seq_along(y), test_idx)
    assert_no_overlap(rownames(x)[train_idx], rownames(x)[test_idx])
    x_train <- x[train_idx, , drop = FALSE]; y_train <- y[train_idx]
    best_cost <- tune_inner(x_train, y_train)
    feats <- fs_select(x_train, y_train, TOP_K)
    std   <- fit_standardizer(x_train[, feats, drop = FALSE])
    x_tr  <- apply_standardizer(x_train[, feats, drop = FALSE], std)
    x_te  <- apply_standardizer(x[test_idx, feats, drop = FALSE], std)
    model <- train_smo_svm(x_tr, y_train, kernel = KERNEL, cost = best_cost,
                           class_weights = TRUE, probability = TRUE)
    pr <- predict_smo_svm(model, x_te, positive_class = POSITIVE)
    preds[[f]] <- data.frame(
      sample_id = rownames(x)[test_idx], truth = as.character(y[test_idx]),
      predicted = pr$predicted_class, prob_pos = pr$probability_positive,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, preds)
}

guarded_summarize <- function(pred_df) {
  m <- compute_binary_metrics(pred_df$truth, pred_df$predicted, pred_df$prob_pos, POSITIVE)
  c(auroc = m$auroc, pr_auc = m$pr_auc, balanced_accuracy = m$balanced_accuracy, mcc = m$mcc,
    sensitivity = sum(pred_df$truth == POSITIVE & pred_df$predicted == POSITIVE) / sum(pred_df$truth == POSITIVE),
    specificity = sum(pred_df$truth != POSITIVE & pred_df$predicted != POSITIVE) / sum(pred_df$truth != POSITIVE))
}

run_guarded <- function(x, y) guarded_summarize(run_nested(x, y))

# ---- one permutation (both arms) -> single result row -----------------------
perm_worker <- function(i, x, y, perm_index) {
  pid <- i - 1L; tstart <- Sys.time()
  y_perm <- y[perm_index[[i]]]
  lk <- run_leaky(x, y_perm); gd <- run_guarded(x, y_perm)
  data.frame(
    permutation_id = pid, is_identity = (pid == 0L),
    leaky_auroc = lk["auroc"], leaky_pr_auc = lk["pr_auc"],
    leaky_balanced_accuracy = lk["balanced_accuracy"], leaky_mcc = lk["mcc"],
    leaky_sensitivity = lk["sensitivity"], leaky_specificity = lk["specificity"],
    guarded_auroc = gd["auroc"], guarded_pr_auc = gd["pr_auc"],
    guarded_balanced_accuracy = gd["balanced_accuracy"], guarded_mcc = gd["mcc"],
    guarded_sensitivity = gd["sensitivity"], guarded_specificity = gd["specificity"],
    gap_auroc = lk["auroc"] - gd["auroc"], gap_pr_auc = lk["pr_auc"] - gd["pr_auc"],
    runtime_secs = as.numeric(difftime(Sys.time(), tstart, units = "secs")),
    stringsAsFactors = FALSE)
}

# ---- driver -----------------------------------------------------------------
main <- function() {
  t0 <- Sys.time()
  message(sprintf("[cfg] RUN_TAG=%s B_PERM=%d FAST=%s PARALLEL=%s N_WORKERS=%d",
                  RUN_TAG, B_PERM, USE_FAST_SELECTOR, USE_PARALLEL, N_WORKERS))
  dd <- load_or_cache(); x <- dd$x; y <- dd$y; n <- length(y)
  message(sprintf("[perm] N=%d; class balance %s", n,
                  paste(names(table(y)), table(y), sep = "=", collapse = ", ")))

  # selector validation (only meaningful when the fast selector is used)
  val <- NULL
  if (isTRUE(USE_FAST_SELECTOR)) {
    message("[validate] comparing original vs fast selector (identity + 5 shuffles)...")
    val <- validate_selector(x, y, n_check = 5)
    print(val)
  }

  # precompute label permutations (serial, reproducible): id 0 = identity
  perm_index <- vector("list", B_PERM + 1L); perm_index[[1]] <- seq_len(n)
  for (b in seq_len(B_PERM)) { set.seed(SEED + PERM_SEED_OFFSET + b); perm_index[[b + 1L]] <- sample.int(n) }

  ids <- seq_along(perm_index)
  run_one <- function(i) perm_worker(i, x, y, perm_index)
  if (isTRUE(USE_PARALLEL) && .Platform$OS.type != "windows" && N_WORKERS > 1) {
    message(sprintf("[run] parallel across permutations (mclapply, %d workers).", N_WORKERS))
    rows <- parallel::mclapply(ids, run_one, mc.cores = N_WORKERS, mc.preschedule = FALSE)
  } else {
    message("[run] serial across permutations.")
    rows <- lapply(ids, run_one)
  }
  nulls <- do.call(rbind, rows)
  nulls <- nulls[order(nulls$permutation_id), ]
  write.csv(nulls, file.path(RESULTS_DIR, sprintf("permutation_%s_null_distributions.csv", RUN_TAG)), row.names = FALSE)

  # identity reproduction check
  idn <- nulls[nulls$permutation_id == 0L, ]
  repro_ok <- abs(idn$leaky_auroc - REF_LEAKY_AUROC) < 0.01 &&
              abs(idn$guarded_auroc - REF_GUARDED_AUROC) < 0.01

  # p-values / null summaries (shuffled only); reference per statistic
  sh <- nulls[nulls$permutation_id > 0L, ]
  prevalence <- mean(y == POSITIVE)
  emp_p <- function(nv, obs) (1 + sum(nv >= obs)) / (length(nv) + 1)
  stat_row <- function(name, nv, obs, ref) data.frame(
    statistic = name, observed = round(obs, 4), reference_value = round(ref, 4),
    null_mean = round(mean(nv), 4), null_sd = round(sd(nv), 4),
    null_p2.5 = round(quantile(nv, 0.025), 4), null_p97.5 = round(quantile(nv, 0.975), 4),
    frac_null_above_reference = round(mean(nv > ref), 4),
    p_obs_ge_null = round(emp_p(nv, obs), 4), n_perm = nrow(sh), stringsAsFactors = FALSE)
  pv <- rbind(
    stat_row("leaky_auroc",    sh$leaky_auroc,    idn$leaky_auroc,    0.5),
    stat_row("guarded_auroc",  sh$guarded_auroc,  idn$guarded_auroc,  0.5),
    stat_row("gap_auroc",      sh$gap_auroc,      idn$gap_auroc,      0),
    stat_row("leaky_pr_auc",   sh$leaky_pr_auc,   idn$leaky_pr_auc,   prevalence),
    stat_row("guarded_pr_auc", sh$guarded_pr_auc, idn$guarded_pr_auc, prevalence),
    stat_row("gap_pr_auc",     sh$gap_pr_auc,     idn$gap_pr_auc,     0))
  write.csv(pv, file.path(RESULTS_DIR, sprintf("permutation_%s_pvalues.csv", RUN_TAG)), row.names = FALSE)

  # figure
  for (dev_fun in list(
    function() png(file.path(FIG_DIR, sprintf("permutation_%s_null.png", RUN_TAG)), width = 1100, height = 500),
    function() pdf(file.path(FIG_DIR, sprintf("permutation_%s_null.pdf", RUN_TAG)), width = 11, height = 5))) {
    dev_fun(); par(mfrow = c(1, 2))
    # Shared breaks/x-range/y-range so BOTH null histograms (leaky ~0.7-1.0 and
    # guarded ~0.5) and the 0.5 reference are fully visible in the same panel.
    auroc_rng <- range(c(sh$leaky_auroc, sh$guarded_auroc, idn$leaky_auroc, idn$guarded_auroc, 0.5))
    pad  <- diff(auroc_rng) * 0.03
    brks <- seq(auroc_rng[1] - pad, auroc_rng[2] + pad, length.out = 26)
    h_leaky   <- hist(sh$leaky_auroc,   breaks = brks, plot = FALSE)
    h_guarded <- hist(sh$guarded_auroc, breaks = brks, plot = FALSE)
    ymax <- max(h_leaky$counts, h_guarded$counts)
    plot(h_leaky, col = "#E99695", main = "Null AUROC under shuffled labels",
         xlab = "AUROC", xlim = range(brks), ylim = c(0, ymax))
    plot(h_guarded, col = "#C2E0C6", add = TRUE)
    abline(v = 0.5, lty = 3); abline(v = idn$leaky_auroc, col = "red", lwd = 2); abline(v = idn$guarded_auroc, col = "darkgreen", lwd = 2)
    legend("topright", bty = "n", cex = 0.8,
           fill = c("#E99695", "#C2E0C6"), legend = c("leaky null", "guarded null"))
    legend("top", bty = "n", cex = 0.8, lwd = 2, col = c("red", "darkgreen"),
           legend = c("observed leaky", "observed guarded"))
    hist(sh$gap_auroc, breaks = 15, col = "#D9E2F3", main = "Null workflow contrast (naive - guarded)", xlab = "delta AUROC")
    abline(v = 0, lty = 3); abline(v = idn$gap_auroc, col = "blue", lwd = 2)
    dev.off()
  }

  # runtime: per-permutation compute time + wall-clock throughput
  mean_compute <- mean(sh$runtime_secs)
  total_wall   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  wall_per_perm <- total_wall / nrow(nulls)
  OLD_PER_PERM <- 361.8
  est <- function(B) sprintf("%.1f min (%.2f h)", B * wall_per_perm / 60, B * wall_per_perm / 3600)

  notes <- c(
    sprintf("# Permutation Control - %s run (GSE25055 only)", RUN_TAG), "",
    sprintf("- Permutations: identity + %d shuffled. Seed %d (label shuffles offset +%d).", B_PERM, SEED, PERM_SEED_OFFSET),
    sprintf("- Fast Welch selector: %s ; parallel: %s (%d workers).", USE_FAST_SELECTOR, USE_PARALLEL, N_WORKERS),
    "- Folds regenerated per permutation (stratified on shuffled labels); FS rerun per permutation.",
    "- Leaky arm: global t-test top-100 before CV. Guarded arm: t-test inside training folds only.",
    "", "## Identity reproduction check",
    sprintf("- leaky AUROC = %.4f (ref %.4f); guarded AUROC = %.4f (ref %.4f). Within 0.01: %s.",
            idn$leaky_auroc, REF_LEAKY_AUROC, idn$guarded_auroc, REF_GUARDED_AUROC, if (repro_ok) "PASS" else "REVIEW"))
  if (!is.null(val)) notes <- c(notes, "", "## Selector validation (original vs fast Welch)",
    sprintf("- %s: overlap %d/%d, exact_order=%s.", val$case, val$overlap, val$top_k, val$exact_order))
  notes <- c(notes, "", "## Null summary (shuffled labels)",
    sprintf("- leaky AUROC null: mean %.4f (%.4f-%.4f); frac > 0.5 = %.3f.", mean(sh$leaky_auroc), min(sh$leaky_auroc), max(sh$leaky_auroc), mean(sh$leaky_auroc > 0.5)),
    sprintf("- guarded AUROC null: mean %.4f (%.4f-%.4f); frac > 0.5 = %.3f.", mean(sh$guarded_auroc), min(sh$guarded_auroc), max(sh$guarded_auroc), mean(sh$guarded_auroc > 0.5)),
    sprintf("- naive-minus-guarded null contrast: mean %.4f (%.4f-%.4f).", mean(sh$gap_auroc), min(sh$gap_auroc), max(sh$gap_auroc)),
    "", "## Runtime",
    sprintf("- mean compute %.1f s/perm; wall %.1f s/perm with %d workers; total %.1f min.", mean_compute, wall_per_perm, N_WORKERS, total_wall / 60),
    sprintf("- old (serial, original selector) was %.1f s/perm.", OLD_PER_PERM),
    sprintf("- estimated B=200: %s ; B=1000: %s (at current wall throughput).", est(200), est(1000)))
  # Footer: smoke runs are code-validation only; full runs (B=200/B=1000) carry
  # the diagnostic interpretation. Detect smoke by RUN_TAG or small B.
  is_smoke <- grepl("smoke", RUN_TAG, ignore.case = TRUE) || B_PERM < 50
  footer <- if (is_smoke) {
    "_Smoke B is small and intended only for code validation._"
  } else {
    "_This permutation run estimates null behavior under shuffled labels. The leaky null distribution should be interpreted as a diagnostic of feature-selection leakage, not as evidence of biological signal. It supports the interpretation that feature-selection leakage can inflate apparent performance._"
  }
  notes <- c(notes, "", footer)
  writeLines(notes, file.path(RESULTS_DIR, sprintf("permutation_%s_notes.md", RUN_TAG)))

  message(sprintf("[perm] DONE (%s). Identity: %s. compute %.1fs/perm; wall %.1fs/perm (%d workers). B=200 ~ %s; B=1000 ~ %s.",
                  RUN_TAG, if (repro_ok) "PASS" else "REVIEW", mean_compute, wall_per_perm, N_WORKERS, est(200), est(1000)))
  if (!repro_ok) message("[perm] WARNING: identity reproduction FAILED (>0.01). Do not trust the null; debug before B=200.")
  invisible(list(nulls = nulls, pvalues = pv, selector_validation = val, repro_ok = repro_ok))
}

if (sys.nframe() == 0) {
  main()
}
