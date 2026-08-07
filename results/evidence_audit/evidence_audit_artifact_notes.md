# Reproducible Omics Evidence Audit — Artifact Notes

This artifact makes the manuscript's Reproducible Omics Evidence Audit reproducible and machine-readable. It performs no new modeling: all values are read from committed outputs.

## Files

- `tables/evidence_audit/evidence_audit_schema.csv`: column dictionary.
- `tables/evidence_audit/evidence_audit_final.csv`: 27 evidence rows across 10 audit domains.
- `tables/evidence_audit/evidence_audit_machine_readable.json`: the same evidence grouped by domain.
- `scripts/10_build_evidence_audit_artifact.R`: canonical generator.

## Canonical sources

- `tables/uncertainty/bootstrap_ci.csv` and `tables/uncertainty/delta_auroc_prauc_ci.csv`
- `results/uncertainty/delong_tests.csv`
- `results/permutation/permutation_b1000_fixed_pvalues.csv`
- `tables/repeated_cv/leakage_gap_by_seed.csv`, `tables/repeated_cv/stability_by_seed.csv`, and `results/repeated_cv/repeated_cv_gap_tests.csv`
- `tables/pilot_gse25055/feature_stability_summary.csv`
- `results/external_validation_gse25065/gse25065_external_metrics.csv`
- `results/external_validation_gse41998/gse41998_metrics.csv` and `tables/external_validation_gse41998/gse41998_bootstrap_ci.csv`

## Interpretation guardrails

The naive and guarded arms differ in feature-selection placement, cost policy, and resampling. Their difference is therefore a descriptive whole-workflow contrast, not an isolated causal estimate of leakage. The paired bootstrap gives ΔAUROC 0.0799 (95% CI 0.0103–0.1552; p=0.028) and ΔPR-AUC 0.0789 (95% CI −0.0426–0.1778; p=0.191). The 30-seed direction check is split-level robustness, not patient-level inference. External results audit transportability and do not validate a clinical biomarker.

Regenerate with `Rscript scripts/10_build_evidence_audit_artifact.R`.
