# External Validation Notes: GSE25055 (discovery) -> GSE25065 (external)

Task 5 of the leakage-aware SMO/SVM omics study. This is the first **external
transportability** test: a model trained and frozen on the discovery cohort
(GSE25055) is applied once to an independent external cohort (GSE25065).

Reproduce with: `source("scripts/05_external_validation_gse25065.R"); main()`
(`R/00_config.R` sets `SEED = 20260620`; fold creation and feature selection are
deterministic under that seed).

## Cohorts

| Role | Accession | Platform | Total | NA excluded | Kept | RD | pCR |
|------|-----------|----------|-------|-------------|------|----|-----|
| Discovery / training | GSE25055 | GPL96 (HG-U133A) | 310 | 4 | 306 | 249 | 57 |
| External validation  | GSE25065 | GPL96 (HG-U133A) | 198 | 16 | 182 | 140 | 42 |

Both cohorts use the label field `pathologic_response_pcr_rd` (pCR vs RD), loaded
by identical logic (same NA handling, same factor levels `c("RD","pCR")`,
positive class = `pCR`). The verified counts reproduce `docs/dataset_audit_report.md`
exactly. GSE25055 and GSE25065 are non-overlapping by design (the original
Hatzis et al. discovery/validation split), and a defensive GSM-id disjointness
check (`assert_no_overlap`) passed.

## What was fit on discovery ONLY, then frozen

The following were estimated using GSE25055 alone and frozen before GSE25065 was
loaded for prediction:

- **Feature selection**: t-test, top K = 100 probes, fit on the full discovery
  cohort (identical strategy to the pilot).
- **SVM cost**: selected by discovery-only 5-fold guarded cross-validation over
  the grid {0.25, 1, 4} (feature selection and scaling re-fit inside each CV
  training split, so the cost choice is not optimistically biased). Selected
  cost = **4** in the submission-aligned, accession-sorted run.
- **Scaling parameters**: center/scale fit on the discovery cohort, selected
  features only.
- **Model**: linear SVM/SMO (`e1071::svm`) with class weights, trained on the
  full discovery cohort.
- **Decision rule**: libsvm default decision (0.5 on calibrated P(pCR)). No
  post-hoc threshold tuning was performed on either cohort.

GSE25065 was **not** used for feature selection, hyperparameter tuning,
scaling-parameter estimation, threshold selection, class balancing, or any
model-design decision. It was scored exactly once.

## Feature alignment

GSE25055 and GSE25065 share platform GPL96, so the frozen probe set transfers
directly: **100 / 100** frozen features were present in GSE25065 (**0 missing**).
No cross-platform mapping or feature imputation was required, and no refitting on
GSE25065 was needed to align. The frozen GSE25055-derived center/scale were
applied as-is.

## External-validation results (GSE25065, n = 182)

| Metric | Value |
|--------|-------|
| AUROC | 0.5633 |
| PR-AUC | 0.2986 |
| Balanced accuracy | 0.5310 |
| MCC | 0.1013 |
| Sensitivity (pCR recall) | 0.1190 |
| Specificity (RD recall) | 0.9429 |

Confusion matrix (rows = truth, cols = predicted):

|        | pred pCR | pred RD |
|--------|----------|---------|
| **pCR** | 5 (TP)  | 37 (FN) |
| **RD**  | 8 (FP)  | 132 (TN)|

## Interpretation

External AUROC (0.563) is below the discovery guarded nested-CV estimate
(0.703 from `results/pilot_gse25055/nested_smo_svm_metrics.csv`), i.e. performance drops
on a truly independent cohort. The frozen model is strongly biased toward the
majority class (RD): specificity is high (0.943) while sensitivity for pCR is low
(0.119), so it identifies residual disease well but misses most pathologic
complete responders. The base rate of pCR in GSE25065 (42 / 182 = 23%) bounds
PR-AUC interpretation; the observed PR-AUC (0.306) is modestly above that base
rate, indicating limited discrimination in this same-study-family cohort.

## Warnings and limitations

- **Modest external performance.** A single-cohort external AUROC near 0.6 is a
  weak transportability result; it should be read as a leakage-honest estimate,
  not as evidence of a deployable biomarker.
- **Low pCR sensitivity.** With the frozen default decision rule the model rarely
  predicts pCR. Threshold re-calibration could trade specificity for sensitivity,
  but any such threshold must be chosen on discovery (or a separate calibration
  set), never on GSE25065, to preserve external validity.
- **External-cohort scope.** GSE25065 is the same-platform, same-study-family audit;
  GSE41998 is reported separately as a cross-platform sensitivity analysis. The
  MDACC sensitivity cohorts remain held out because of overlap risk.
- **Same platform / same study family.** GSE25055 and GSE25065 are from the same
  authors and platform; this result alone does not establish broader transportability.
- **No probe-to-gene collapsing.** Analysis is at the probe level on GPL96; a
  gene-level summarization could change which features transfer.
- **Class imbalance handled only via SVM class weights**; SMOTE and alternative
  imbalance strategies were not used (consistent with the pilot).
- Results are deterministic under `SEED = 20260620`; different seeds change CV
  fold assignment and could shift the selected cost and downstream metrics.
