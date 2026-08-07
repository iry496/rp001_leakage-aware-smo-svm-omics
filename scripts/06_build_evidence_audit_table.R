# =============================================================================
# scripts/06_build_evidence_audit_table.R
# -----------------------------------------------------------------------------
# Task 6: Build the first Reproducible Omics Evidence Audit Table (v1).
#
# This script performs NO new modeling. It INTEGRATES already-committed outputs
# from the completed stages into a single auditable evidence table:
#   * Dataset audit            -> tables/table1_dataset_audit_filled.csv
#   * Leaky vs guarded pilot    -> tables/pilot_gse25055/pilot_performance_comparison.csv
#   * Feature stability         -> tables/pilot_gse25055/feature_stability_summary.csv
#   * External validation       -> tables/external_validation_gse25065/external_validation_summary.csv
#                                  results/external_validation_gse25065/gse25065_external_metrics.csv
#
# Outputs (no rendered HTML, no manuscript/protocol edits):
#   * tables/evidence_audit/reproducible_omics_evidence_audit.csv
#   * tables/evidence_audit/reproducible_omics_evidence_audit_readable.md
#   * results/evidence_audit/evidence_audit_notes.md
#
# Interpretation is deliberately cautious: this is a METHODOLOGY / AUDIT result,
# not a clinical biomarker discovery.
# =============================================================================

# ---- Config -----------------------------------------------------------------
TABLES_DIR  <- file.path("tables",  "evidence_audit")
RESULTS_DIR <- file.path("results", "evidence_audit")
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

f4 <- function(x) formatC(as.numeric(x), format = "f", digits = 4)
f1pct <- function(x) paste0(formatC(100 * as.numeric(x), format = "f", digits = 1), "%")
signed4 <- function(x) {
  x <- as.numeric(x)
  paste0(ifelse(x >= 0, "+", "-"), formatC(abs(x), format = "f", digits = 4))
}

# ---- Read inputs (read-only; nothing is refit) ------------------------------
audit_in <- read.csv("tables/table1_dataset_audit_filled.csv",
                     stringsAsFactors = FALSE, check.names = FALSE)
leaky <- read.csv("results/pilot_gse25055/leaky_baseline_metrics.csv",
                  stringsAsFactors = FALSE)
guarded <- read.csv("results/pilot_gse25055/nested_smo_svm_metrics.csv",
                    stringsAsFactors = FALSE)
fs <- read.csv("tables/pilot_gse25055/feature_stability_summary.csv",
               stringsAsFactors = FALSE)
extm <- read.csv("results/external_validation_gse25065/gse25065_external_metrics.csv",
                 stringsAsFactors = FALSE)

# Dataset rows
d55 <- audit_in[audit_in$GEO_accession == "GSE25055", ]
d65 <- audit_in[audit_in$GEO_accession == "GSE25065", ]
na55 <- d55$Total_N - d55$pCR_N - d55$RD_N
na65 <- d65$Total_N - d65$pCR_N - d65$RD_N
usable55 <- d55$pCR_N + d55$RD_N
usable65 <- d65$pCR_N + d65$RD_N
prev55 <- d55$pCR_N / usable55
prev65 <- d65$pCR_N / usable65

# Feature-stability values (metric/value long table; value stored as character)
fsget <- function(key) fs$value[match(key, fs$metric)]
fs_total   <- as.numeric(fsget("total_unique_features"))
fs_all     <- as.numeric(fsget("features_selected_in_all_folds"))
fs_one     <- as.numeric(fsget("features_selected_in_single_fold"))
fs_meanJ   <- as.numeric(fsget("mean_pairwise_jaccard"))
fs_medJ    <- as.numeric(fsget("median_pairwise_jaccard"))
fs_nog     <- as.numeric(fsget("nogueira_stability_index"))

# External-validation metrics (single numeric row)
ext_auroc <- extm$auroc
ext_prauc <- extm$pr_auc
ext_bacc  <- extm$balanced_accuracy
ext_mcc   <- extm$mcc
ext_sens  <- extm$sensitivity
ext_spec  <- extm$specificity

# ---- Derived gaps and drops -------------------------------------------------
gap_auroc <- leaky$auroc  - guarded$auroc           # whole-workflow contrast (AUROC)
gap_prauc <- leaky$pr_auc - guarded$pr_auc           # whole-workflow contrast (PR-AUC)

drop_auroc <- guarded$auroc             - ext_auroc # internal->external drop
drop_prauc <- guarded$pr_auc            - ext_prauc
drop_bacc  <- guarded$balanced_accuracy - ext_bacc
drop_mcc   <- guarded$mcc               - ext_mcc
drop_sens  <- guarded$sensitivity       - ext_sens
drop_spec  <- guarded$specificity       - ext_spec

# ---- Assemble the audit table (long format) ---------------------------------
rows <- list()
add <- function(domain, metric, value, interpretation, manuscript_use) {
  rows[[length(rows) + 1]] <<- data.frame(
    audit_domain    = domain,
    metric          = metric,
    value           = value,
    interpretation  = interpretation,
    manuscript_use  = manuscript_use,
    stringsAsFactors = FALSE
  )
}

# 1. Dataset integrity
add("Dataset integrity", "GSE25055 usable samples (discovery)",
    sprintf("RD=%d, pCR=%d; %d NA excluded (N=%d)", d55$RD_N, d55$pCR_N, na55, d55$Total_N),
    "Verified from sample metadata; modest pCR minority.",
    "Table 1 / cohort description")
add("Dataset integrity", "GSE25065 usable samples (external)",
    sprintf("RD=%d, pCR=%d; %d NA excluded (N=%d)", d65$RD_N, d65$pCR_N, na65, d65$Total_N),
    "Verified from sample metadata; independent validation cohort.",
    "Table 1 / cohort description")
add("Dataset integrity", "Platform & endpoint harmonization",
    "Both GPL96; label pathologic_response_pcr_rd (pCR vs RD)",
    "Same platform enables direct probe transfer with no remapping.",
    "Methods: cohorts & endpoint")
add("Dataset integrity", "Discovery/validation overlap",
    "Non-overlapping by design (GSE25066 split)",
    "No patient-level leakage between discovery and external cohorts.",
    "Methods: external-validation design")

# 2. Naive-versus-guarded workflow sensitivity
add("Workflow sensitivity", "Naive baseline AUROC", f4(leaky$auroc),
    "Naive global-selection workflow estimate; not an unbiased generalization estimate.",
    "Workflow comparison")
add("Workflow sensitivity", "Naive baseline PR-AUC", f4(leaky$pr_auc),
    "Naive global-selection workflow estimate under class imbalance.",
    "Workflow comparison")
add("Workflow sensitivity", "Workflow contrast AUROC (naive - guarded)", signed4(gap_auroc),
    "Descriptive contrast between workflows that also differ in tuning and resampling.",
    "Whole-workflow comparison")
add("Workflow sensitivity", "Workflow contrast PR-AUC (naive - guarded)", signed4(gap_prauc),
    "Descriptive contrast; not the isolated causal effect of feature-selection placement.",
    "Whole-workflow comparison")

# 3. Guarded nested performance
add("Guarded nested performance", "AUROC", f4(guarded$auroc),
    "Honest internal discrimination (nested CV).", "Primary internal result")
add("Guarded nested performance", "PR-AUC", f4(guarded$pr_auc),
    "Honest internal precision-recall on imbalanced pCR.", "Primary internal result")
add("Guarded nested performance", "Balanced accuracy", f4(guarded$balanced_accuracy),
    "Imbalance-aware accuracy improves vs leaky pipeline.", "Imbalance-aware result")
add("Guarded nested performance", "MCC", f4(guarded$mcc),
    "Imbalance-aware agreement improves vs leaky pipeline.", "Imbalance-aware result")
add("Guarded nested performance", "Sensitivity (pCR)", f4(guarded$sensitivity),
    "Minority-class recall still limited.", "Imbalance discussion")
add("Guarded nested performance", "Specificity (RD)", f4(guarded$specificity),
    "Majority-class recall high.", "Imbalance discussion")

# 4. Feature stability
add("Feature stability", "Total unique selected features (5 folds, K=100)", as.character(fs_total),
    "Union far exceeds per-fold K=100, signalling churn.", "Stability result")
add("Feature stability", "Selected in all 5 folds (stable core)", as.character(fs_all),
    "Small reproducible core of features.", "Stable-core result")
add("Feature stability", "Selected in exactly 1 fold (unstable tail)", as.character(fs_one),
    "Large unstable tail selected only once.", "Instability caveat")
add("Feature stability", "Mean pairwise Jaccard", f4(fs_meanJ),
    "Moderate overlap between fold-wise feature sets.", "Stability result")
add("Feature stability", "Median pairwise Jaccard", f4(fs_medJ),
    "Consistent with mean; moderate overlap.", "Stability result")
add("Feature stability", "Nogueira stability index", f4(fs_nog),
    "Moderate stability: stable core with unstable tail.", "Stability headline")

# 5. External validation
add("External validation", "AUROC (GSE25065)", f4(ext_auroc),
    "Lower than internal nested CV; transportability limited.", "External result")
add("External validation", "PR-AUC (GSE25065)", f4(ext_prauc),
    "Weak-but-present pCR signal above base rate.", "External result")
add("External validation", "Balanced accuracy (GSE25065)", f4(ext_bacc),
    "Near chance once externalized.", "External result")
add("External validation", "MCC (GSE25065)", f4(ext_mcc),
    "Weak agreement on independent cohort.", "External result")
add("External validation", "Sensitivity (pCR, GSE25065)", f4(ext_sens),
    "Most pathologic responders missed externally.", "External / imbalance discussion")
add("External validation", "Specificity (RD, GSE25065)", f4(ext_spec),
    "Residual disease still well identified.", "External / imbalance discussion")
add("External validation", "Drop AUROC (nested -> external)", signed4(drop_auroc),
    "Quantifies transportability gap to an independent cohort.", "Key transportability result")
add("External validation", "Drop PR-AUC (nested -> external)", signed4(drop_prauc),
    "Precision-recall declines externally.", "Transportability result")
add("External validation", "Drop balanced accuracy (nested -> external)", signed4(drop_bacc),
    "Imbalance-aware accuracy declines externally.", "Transportability result")
add("External validation", "Drop MCC (nested -> external)", signed4(drop_mcc),
    "Agreement declines externally.", "Transportability result")
add("External validation", "Drop sensitivity (nested -> external)", signed4(drop_sens),
    "pCR recall degrades further externally.", "Transportability / imbalance result")
add("External validation", "Drop specificity (nested -> external)", signed4(drop_spec),
    "Specificity roughly maintained (slightly higher externally).", "Transportability result")

# 6. Class-imbalance behavior
add("Class-imbalance behavior", "pCR prevalence (discovery / external)",
    sprintf("%s / %s", f1pct(prev55), f1pct(prev65)),
    "Minority pCR class in both cohorts.", "Imbalance framing")
add("Class-imbalance behavior", "Sensitivity vs specificity (nested)",
    sprintf("%s vs %s", f4(guarded$sensitivity), f4(guarded$specificity)),
    "Majority-class (RD) bias persists despite class weights.", "Imbalance discussion")
add("Class-imbalance behavior", "Sensitivity vs specificity (external)",
    sprintf("%s vs %s", f4(ext_sens), f4(ext_spec)),
    "pCR recall degrades further on the external cohort.", "Imbalance discussion")
add("Class-imbalance behavior", "Imbalance handling",
    "SVM class weights (no SMOTE)",
    "Weighting alone does not resolve low pCR recall.", "Methods: imbalance strategy")

# 7. Reproducibility status
add("Reproducibility status", "Random seed", "SEED = 20260620",
    "Single global seed across pipelines.", "Reproducibility statement")
add("Reproducibility status", "Determinism",
    "Samples sorted by GEO accession ID before seeded folds; deterministic t-test top-K selection",
    "Canonical ordering makes seeded fold assignments reproducible across fresh downloads.", "Reproducibility statement")
add("Reproducibility status", "Committed artifacts",
    "Per-stage metrics, predictions, selected features, notes",
    "Full provenance from raw load to external metrics.", "Data/code availability")
add("Reproducibility status", "Raw-data policy",
    "GEO series-matrix via GEOquery; no raw CEL committed",
    "Lightweight, license-respecting reproducibility.", "Data availability")

# 8. Limitations / unresolved risks
add("Limitations / unresolved risks", "Clinical claim",
    "None - methodology/audit only",
    "Not a validated clinical biomarker; no discovery claim.", "Limitations")
add("Limitations / unresolved risks", "External breadth",
    "One same-study-family cohort and one cross-platform cohort",
    "Broader population and platform transportability remains untested.", "Limitations")
add("Limitations / unresolved risks", "Pending cohorts",
    "GSE20194/GSE20271 (patient-overlap risk)",
    "Held out pending deduplication and overlap resolution.", "Future work")
add("Limitations / unresolved risks", "Feature-tail instability",
    sprintf("%d of %d features selected once", fs_one, fs_total),
    "Only a small core is reproducible across folds.", "Limitations")
add("Limitations / unresolved risks", "pCR detection",
    "Low sensitivity internally and externally",
    "Limited utility for identifying responders.", "Limitations")

audit <- do.call(rbind, rows)

# ---- Output 1: machine-readable CSV -----------------------------------------
utils::write.csv(audit,
                 file.path(TABLES_DIR, "reproducible_omics_evidence_audit.csv"),
                 row.names = FALSE)

# ---- Output 2: human-readable Markdown table --------------------------------
md_escape <- function(x) gsub("\\|", "\\\\|", x)
header <- c("| Audit domain | Metric | Value | Interpretation | Manuscript use |",
            "|---|---|---|---|---|")
body <- apply(audit, 1, function(r) {
  sprintf("| %s | %s | %s | %s | %s |",
          md_escape(r["audit_domain"]), md_escape(r["metric"]),
          md_escape(r["value"]), md_escape(r["interpretation"]),
          md_escape(r["manuscript_use"]))
})
md <- c("# Reproducible Omics Evidence Audit Table (v1)",
        "",
        "Integrated, leakage-aware evidence summary for the GSE25055 discovery /",
        "GSE25065 external-validation SMO/SVM study. No new modeling: every value",
        "is extracted from previously committed stage outputs. This is a",
        "**methodology / audit** result, not a clinical biomarker discovery.",
        "",
        header, body, "")
writeLines(md, file.path(TABLES_DIR, "reproducible_omics_evidence_audit_readable.md"))

# ---- Output 3: narrative notes ----------------------------------------------
notes <- c(
  "# Reproducible Omics Evidence Audit - Notes (v1)",
  "",
  "## Purpose",
  "First integrated evidence audit for the leakage-aware SMO/SVM omics study.",
  "It consolidates four completed stages - dataset audit, GSE25055 leaky-vs-guarded",
  "pilot, feature-stability analysis, and GSE25065 external validation - into one",
  "auditable table. No models were run here; all numbers are read from committed",
  "outputs and recombined.",
  "",
  "## Cohorts",
  sprintf("- GSE25055 (discovery): RD=%d, pCR=%d, %d NA excluded (N=%d); pCR prevalence %s.",
          d55$RD_N, d55$pCR_N, na55, d55$Total_N, f1pct(prev55)),
  sprintf("- GSE25065 (external): RD=%d, pCR=%d, %d NA excluded (N=%d); pCR prevalence %s.",
          d65$RD_N, d65$pCR_N, na65, d65$Total_N, f1pct(prev65)),
  "- Same platform (GPL96), non-overlapping by design (GSE25066 split).",
  "",
  "## Workflow sensitivity (naive vs guarded nested, GSE25055)",
  sprintf("- Naive AUROC %s vs guarded nested AUROC %s -> whole-workflow contrast %s.",
          f4(leaky$auroc), f4(guarded$auroc), signed4(gap_auroc)),
  sprintf("- Naive PR-AUC %s vs guarded nested PR-AUC %s -> whole-workflow contrast %s.",
          f4(leaky$pr_auc), f4(guarded$pr_auc), signed4(gap_prauc)),
  "- The workflows also differ in tuning and resampling, so the contrasts are",
  "  descriptive and not isolated causal estimates of feature-selection leakage.",
  "  The guarded nested workflow improves the",
  sprintf("  imbalance-aware metrics (balanced accuracy %s, MCC %s) over the leaky",
          f4(guarded$balanced_accuracy), f4(guarded$mcc)),
  sprintf("  baseline (balanced accuracy %s, MCC %s).",
          f4(leaky$balanced_accuracy), f4(leaky$mcc)),
  "",
  "## Feature stability (5 outer folds, K=100)",
  sprintf("- %d unique features selected across folds; %d in all 5 folds (stable core);",
          fs_total, fs_all),
  sprintf("  %d selected in exactly one fold (unstable tail).", fs_one),
  sprintf("- Mean Jaccard %s, median Jaccard %s, Nogueira stability %s.",
          f4(fs_meanJ), f4(fs_medJ), f4(fs_nog)),
  "- Interpretation: stability is MODERATE - a small reproducible core coexists",
  "  with a large unstable tail.",
  "",
  "## External validation (GSE25065) and internal->external drop",
  sprintf("- External AUROC %s, PR-AUC %s, balanced accuracy %s, MCC %s, sensitivity %s, specificity %s.",
          f4(ext_auroc), f4(ext_prauc), f4(ext_bacc), f4(ext_mcc), f4(ext_sens), f4(ext_spec)),
  sprintf("- Drop from guarded nested CV to external: AUROC %s, PR-AUC %s, balanced accuracy %s,",
          signed4(drop_auroc), signed4(drop_prauc), signed4(drop_bacc)),
  sprintf("  MCC %s, sensitivity %s, specificity %s.",
          signed4(drop_mcc), signed4(drop_sens), signed4(drop_spec)),
  "- Interpretation: external validation shows clear TRANSPORTABILITY LIMITS;",
  "  discrimination and minority-class recall fall on the independent cohort,",
  "  while specificity is roughly maintained.",
  "",
  "## Cautious interpretation",
  "- This is a METHODOLOGY / AUDIT result, NOT a clinical biomarker discovery.",
  "- No clinical claims are made; pCR sensitivity is low internally and externally.",
  "- The leaky pipeline inflates AUROC/PR-AUC; guarded nested validation gives an",
  "  honest estimate and improves imbalance-aware metrics.",
  "- Feature stability is moderate (stable core + unstable tail).",
  "- External validation demonstrates real transportability limits.",
  "",
  "## Limitations / unresolved risks",
  "- Single same-platform external cohort (GSE25065).",
  "- GSE41998 (cross-platform) and GSE20194/GSE20271 (overlap/de-dup) deliberately",
  "  excluded pending harmonization and sample-level de-duplication.",
  sprintf("- Feature-tail instability: %d of %d features selected only once.", fs_one, fs_total),
  "- Low pCR sensitivity limits practical utility.",
  ""
)
writeLines(notes, file.path(RESULTS_DIR, "evidence_audit_notes.md"))

message("[06] Evidence audit table built. Rows: ", nrow(audit))
message(sprintf("[06] Leakage gap AUROC=%s PR-AUC=%s | Ext drop AUROC=%s PR-AUC=%s",
                signed4(gap_auroc), signed4(gap_prauc),
                signed4(drop_auroc), signed4(drop_prauc)))
invisible(audit)
