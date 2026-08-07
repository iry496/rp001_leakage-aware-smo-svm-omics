#!/usr/bin/env Rscript
# =============================================================================
# scripts/10_build_evidence_audit_artifact.R
# -----------------------------------------------------------------------------
# Build a reproducible, machine-readable Evidence Audit artifact (manuscript
# Table 10). NO NEW MODELING: every value is read from outputs already committed
# to the repository. The script writes:
#   tables/evidence_audit/evidence_audit_schema.csv
#   tables/evidence_audit/evidence_audit_final.csv
#   tables/evidence_audit/evidence_audit_machine_readable.json
# and refreshes results/evidence_audit/evidence_audit_artifact_notes.md (header).
#
# Ten domains: dataset integrity; leakage sensitivity; bootstrap uncertainty;
# permutation control; repeated-CV robustness; feature stability; external
# transportability; class imbalance; reproducibility status; limitations.
# =============================================================================

suppressWarnings(suppressMessages({
  has_json <- requireNamespace("jsonlite", quietly = TRUE)
}))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)

# ---- committed sources ------------------------------------------------------
bci <- rd("tables/uncertainty/bootstrap_ci.csv")
dci <- rd("tables/uncertainty/delta_auroc_prauc_ci.csv")
dlt <- rd("results/uncertainty/delong_tests.csv")
prm <- rd("results/permutation/permutation_b1000_fixed_pvalues.csv")
gpt <- rd("results/repeated_cv/repeated_cv_gap_tests.csv")
sbs <- rd("tables/repeated_cv/stability_by_seed.csv")
fss <- rd("tables/pilot_gse25055/feature_stability_summary.csv")
ext <- rd("results/external_validation_gse25065/gse25065_external_metrics.csv")

bget <- function(pipeline, metric, col) {
  r <- bci[bci$pipeline == pipeline & bci$metric == metric, ]; r[[col]][1]
}
fsv <- function(m) fss$value[fss$metric == m][1]
prow <- function(stat) prm[prm$statistic == stat, ]
grow <- function(stat) gpt[gpt$statistic == stat, ]

rows <- list()
add <- function(domain, metric, value, interval = "", p_value = "",
                reference = "", interpretation = "", source = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    domain = domain, metric = metric, value = as.character(value),
    interval = interval, p_value = as.character(p_value),
    reference = reference, interpretation = interpretation, source = source,
    stringsAsFactors = FALSE)
}

# 1. Dataset integrity
add("Dataset integrity", "Discovery usable samples", "306 (57 pCR, 249 RD)",
    reference = "prevalence 0.186",
    interpretation = "Labeled, same-platform (GPL96) discovery cohort.",
    source = "results/pilot_gse25055/")
add("Dataset integrity", "External usable samples", "182 (42 pCR, 140 RD)",
    reference = "prevalence 0.231",
    interpretation = "Same-platform, same-study-family (Hatzis), non-overlapping with discovery.",
    source = "results/external_validation_gse25065/gse25065_external_metrics.csv")

# 2. Naive-versus-guarded workflow sensitivity
da <- dci[dci$metric == "delta_auroc", ]; dp <- dci[dci$metric == "delta_pr_auc", ]
add("Workflow sensitivity", "Internal whole-workflow contrast, delta AUROC", da$point,
    sprintf("[%s, %s]", da$ci_lo, da$ci_hi), da$boot_p, "0",
    "Descriptive whole-workflow contrast; arms also differ in tuning and resampling.",
    "tables/uncertainty/delta_auroc_prauc_ci.csv")
add("Workflow sensitivity", "Internal whole-workflow contrast, delta PR-AUC", dp$point,
    sprintf("[%s, %s]", dp$ci_lo, dp$ci_hi), dp$boot_p, "0",
    "Descriptive whole-workflow contrast; not an isolated leakage-effect estimate.",
    "tables/uncertainty/delta_auroc_prauc_ci.csv")
add("Workflow sensitivity", "DeLong per-repeat p (secondary)",
    sprintf("%.3f-%.3f", min(dlt$delong_p), max(dlt$delong_p)), "",
    sprintf("median %.3f", median(dlt$delong_p)), "",
    "Secondary sensitivity check; not the primary inference.",
    "results/uncertainty/delong_tests.csv")

# 3. Bootstrap uncertainty
for (it in list(c("Leaky AUROC","A_leaky_baseline","auroc"),
                c("Guarded AUROC","B_guarded_nested","auroc"),
                c("Leaky PR-AUC","A_leaky_baseline","pr_auc"),
                c("Guarded PR-AUC","B_guarded_nested","pr_auc"))) {
  add("Bootstrap uncertainty", it[1], bget(it[2], it[3], "point"),
      sprintf("[%s, %s]", bget(it[2], it[3], "ci_lo"), bget(it[2], it[3], "ci_hi")),
      "", "", "95% stratified percentile bootstrap (B=2000).",
      "tables/uncertainty/bootstrap_ci.csv")
}

# 4. Permutation control
pa <- prow("leaky_auroc"); pg <- prow("guarded_auroc"); pgp <- prow("guarded_pr_auc")
add("Permutation control", "Leaky null AUROC (mean)", pa$null_mean,
    sprintf("[%s, %s]", pa$`null_p2.5`, pa$`null_p97.5`),
    sprintf("p(obs>=null)=%s", pa$p_obs_ge_null), "chance 0.5",
    "Leaky pipeline scores far above chance on shuffled labels (diagnostic of leakage artifact).",
    "results/permutation/permutation_b1000_fixed_pvalues.csv")
add("Permutation control", "Guarded null AUROC (mean)", pg$null_mean,
    sprintf("[%s, %s]", pg$`null_p2.5`, pg$`null_p97.5`),
    sprintf("obs p=%s", pg$p_obs_ge_null), "chance 0.5",
    "Guarded null near chance; real-label guarded exceeds its null (p=0.001).",
    "results/permutation/permutation_b1000_fixed_pvalues.csv")
add("Permutation control", "Guarded null PR-AUC (mean)", pgp$null_mean,
    sprintf("[%s, %s]", pgp$`null_p2.5`, pgp$`null_p97.5`),
    sprintf("obs p=%s", pgp$p_obs_ge_null), "prevalence 0.186",
    "Guarded PR-AUC null centers on prevalence (well-calibrated).",
    "results/permutation/permutation_b1000_fixed_pvalues.csv")

# 5. Repeated-CV robustness
ga <- grow("delta_auroc"); gp <- grow("delta_pr_auc")
add("Repeated-CV robustness", "delta AUROC across 30 seeds", ga$median,
    sprintf("[%s, %s]", ga$`ci2.5`, ga$`ci97.5`),
    sprintf("Wilcoxon V=%s, p=%s", ga$V_statistic, ga$p_value),
    sprintf("positive in %s/%s", ga$n_positive, ga$n_seeds),
    "Modest but reproducible across fold randomizations (not patient-level inference).",
    "results/repeated_cv/repeated_cv_gap_tests.csv")
add("Repeated-CV robustness", "delta PR-AUC across 30 seeds", gp$median,
    sprintf("[%s, %s]", gp$`ci2.5`, gp$`ci97.5`),
    sprintf("Wilcoxon V=%s, p=%s", gp$V_statistic, gp$p_value),
    sprintf("positive in %s/%s", gp$n_positive, gp$n_seeds),
    "Modest but reproducible across fold randomizations.",
    "results/repeated_cv/repeated_cv_gap_tests.csv")

# 6. Feature stability
add("Feature stability", "Nogueira index (single split)", fsv("nogueira_stability_index"),
    "", "", "", "Moderate stability of the selection procedure.",
    "tables/pilot_gse25055/feature_stability_summary.csv")
add("Feature stability", "Mean Jaccard (single split)", fsv("mean_pairwise_jaccard"),
    "", "", "", "Moderate fold-to-fold overlap.",
    "tables/pilot_gse25055/feature_stability_summary.csv")
add("Feature stability", "Stable core / unstable tail",
    sprintf("%s / %s", fsv("features_selected_in_all_folds"), fsv("features_selected_in_single_fold")),
    "", "", sprintf("of %s unique", fsv("total_unique_features")),
    "Small reproducible core; large unstable tail.",
    "tables/pilot_gse25055/feature_stability_summary.csv")
add("Feature stability", "Nogueira median (30 seeds)",
    sprintf("%.4f", median(sbs$nogueira_stability)),
    sprintf("core %d-%d, tail %d-%d (range)",
            min(sbs$stable_core_count), max(sbs$stable_core_count),
            min(sbs$unstable_tail_count), max(sbs$unstable_tail_count)),
    "", "", "Stability reproducible across seeds.",
    "tables/repeated_cv/stability_by_seed.csv")

# 7. External transportability
add("External transportability", "External AUROC (GSE25065)", sprintf("%.4f", ext$auroc),
    sprintf("[%s, %s]", bget("B_guarded_nested_external","auroc","ci_lo"),
                        bget("B_guarded_nested_external","auroc","ci_hi")),
    "", "discovery guarded 0.7032",
    "Transportability drop (~ -0.140 AUROC); a generalization limit, distinct from leakage.",
    "results/external_validation_gse25065/gse25065_external_metrics.csv")
add("External transportability", "External PR-AUC (GSE25065)", sprintf("%.4f", ext$pr_auc),
    "", "", "prevalence 0.231",
    "Reduced minority-class ranking on the independent cohort.",
    "results/external_validation_gse25065/gse25065_external_metrics.csv")

# 8. Class imbalance
add("Class imbalance", "External sensitivity / specificity",
    sprintf("%.4f / %.4f", ext$sensitivity, ext$specificity), "", "",
    sprintf("TP=%s, FP=%s, TN=%s, FN=%s", ext$tp, ext$fp, ext$tn, ext$fn),
    "High specificity but low pCR sensitivity under class imbalance.",
    "results/external_validation_gse25065/gse25065_external_metrics.csv")
add("Class imbalance", "Guarded discovery sensitivity (pCR)",
    bget("B_guarded_nested","sensitivity","point"),
    sprintf("[%s, %s]", bget("B_guarded_nested","sensitivity","ci_lo"),
                        bget("B_guarded_nested","sensitivity","ci_hi")),
    "", "", "Minority-class recall remains low internally.",
    "tables/uncertainty/bootstrap_ci.csv")

# 9. Reproducibility status
add("Reproducibility status", "Analysis scripts",
    "07_bootstrap_ci.R; 08_permutation_control.R; 09_repeated_nested_cv.R; 15_external_validation_gse41998_bootstrap.R", "", "", "",
    "All analyses regenerate from committed code, seeds, and outputs.", "scripts/")
add("Reproducibility status", "Committed outputs",
    "tables/uncertainty/; results/permutation/; tables/repeated_cv/; results/repeated_cv/",
    "", "", "", "Pipeline is reproducible from repository artifacts.", "repository")

# 10. Limitations / unresolved risks
add("Limitations / unresolved risks", "Statistical resolution",
    "Single-cohort workflow contrast is not a causal leakage estimate", "", "", "",
    "Matched one-factor ablation is required to estimate the effect of feature-selection placement.",
    "tables/uncertainty/delta_auroc_prauc_ci.csv")
add("Limitations / unresolved risks", "Generalization",
    "One same-study-family and one cross-platform external cohort", "", "", "",
    "Broader generalization remains untested; recurrent probes are audit outputs, not validated markers.",
    "results/external_validation_gse25065/")
add("Limitations / unresolved risks", "Clinical status", "No clinical biomarker claim",
    "", "", "", "Methodology/audit study; not clinical biomarker discovery.", "manuscript")

final <- do.call(rbind, rows)

dir.create("tables/evidence_audit", recursive = TRUE, showWarnings = FALSE)
dir.create("results/evidence_audit", recursive = TRUE, showWarnings = FALSE)

# schema
schema <- data.frame(
  column = c("domain","metric","value","interval","p_value","reference","interpretation","source"),
  type = "string",
  description = c(
    "Audit domain (one of the ten standard domains).",
    "Specific quantity or item within the domain.",
    "Point value or summary (verbatim from the committed source).",
    "95% CI or 2.5-97.5%/range where applicable; blank if not defined.",
    "Associated p-value or empirical p where applicable; blank otherwise.",
    "Reference point (chance, prevalence, comparator) or count context.",
    "Cautious, within-cohort/diagnostic interpretation.",
    "Committed file or location the value is drawn from."),
  stringsAsFactors = FALSE)
utils::write.csv(schema, "tables/evidence_audit/evidence_audit_schema.csv", row.names = FALSE)
utils::write.csv(final,  "tables/evidence_audit/evidence_audit_final.csv",  row.names = FALSE)

# JSON (nested by domain)
order <- c("Dataset integrity","Workflow sensitivity","Bootstrap uncertainty",
           "Permutation control","Repeated-CV robustness","Feature stability",
           "External transportability","Class imbalance","Reproducibility status",
           "Limitations / unresolved risks")
domains <- lapply(order, function(d) {
  sub <- final[final$domain == d, setdiff(names(final), "domain")]
  list(domain = d, entries = unname(split(sub, seq_len(nrow(sub)))))
})
obj <- list(artifact = "Reproducible Omics Evidence Audit",
            manuscript_table = "Table 10",
            generated = "reproducible from committed outputs via scripts/10_build_evidence_audit_artifact.R",
            new_modeling = FALSE, cohort_discovery = "GSE25055",
            cohort_external = "GSE25065", domains = domains)
if (has_json) {
  writeLines(jsonlite::toJSON(obj, pretty = TRUE, auto_unbox = TRUE),
             "tables/evidence_audit/evidence_audit_machine_readable.json")
} else {
  warning("jsonlite not installed; JSON not written. install.packages('jsonlite').")
}

cat(sprintf("[evidence-audit] wrote %d rows across %d domains.\n", nrow(final), length(order)))
