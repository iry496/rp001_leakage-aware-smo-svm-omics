# Uncertainty Quantification — Notes

These analyses resample committed out-of-fold predictions; they do not refit models.

## Inputs

- `results/pilot_gse25055/leaky_baseline_predictions.csv`
- `results/pilot_gse25055/nested_smo_svm_predictions.csv`
- `results/external_validation_gse25065/gse25065_external_predictions.csv`

## Method

Stratified percentile bootstrap, B=2000, seed 20260620. The naive-versus-guarded comparison uses a paired patient-level bootstrap. DeLong tests are run per repeat because the naive arm has repeated-CV scores whereas the guarded arm has one out-of-fold score per sample.

## Submission-aligned findings

- Naive AUROC 0.7830; guarded AUROC 0.7032; ΔAUROC 0.0799 (95% CI 0.0103–0.1552; p=0.028).
- Naive PR-AUC 0.4257; guarded PR-AUC 0.3476; ΔPR-AUC 0.0789 (95% CI −0.0426–0.1778; p=0.191).
- The AUROC contrast is resolved as positive in this resampling analysis; the PR-AUC contrast is not.
- The arms also differ in tuning and resampling, so neither contrast isolates the causal effect of feature-selection placement.

The committed CSVs are reproduced by `scripts/checks/07_bootstrap_ci_crosscheck.py`; the canonical R implementation is `scripts/07_bootstrap_ci.R`.
