# Protocol v0.1

> Historical protocol retained for provenance. It predates the final dataset
> roles, manuscript title, and Journal of Biomedical Informatics submission;
> the v1.2.1 manuscript and release supersede its provisional decisions.

## Manuscript title

A Leakage-Aware Data Mining Framework for Transcriptomic Biomarker Classification Using Nested SMO/SVM and Feature-Stability Auditing

## Target journal

BioData Mining, Methodology article.

## Study objective

Quantify how feature-selection leakage affects SMO/SVM transcriptomic biomarker classification and evaluate whether nested validation plus feature-stability auditing produces more realistic, reproducible evidence.

## Primary research questions

1. How much does global feature selection inflate apparent SMO/SVM performance?
2. What is the realistic performance under guarded nested cross-validation?
3. Are selected transcriptomic features stable across folds and resampling?
4. Does the guarded model transport to independent public GEO cohorts?
5. Can a compact evidence-audit table summarize leakage sensitivity, stability, and validation evidence?

## Dataset roles to lock after audit

| Role | Candidate datasets | Decision status |
|---|---|---|
| Discovery / nested CV | GSE25055 or GSE25065 | pending |
| External validation 1 | Non-overlapping Hatzis subseries | pending |
| External validation 2 | GSE20194 | pending |
| Sensitivity | GSE41998 or GSE20271 | pending |

## Endpoint

Binary endpoint: pathological complete response (pCR) versus residual disease (RD). Ambiguous or incompatible response labels must be excluded before modeling.

## Pipeline A: naive / leaky baseline

Purpose: simulate a common flawed workflow.

1. Use the full discovery cohort.
2. Perform supervised feature selection globally using all labels.
3. Select top K features.
4. Cross-validate SMO/SVM on the preselected feature set.
5. Report apparent performance.

## Pipeline B: guarded nested SMO/SVM

Purpose: estimate realistic generalization.

1. Outer loop holds out a test fold.
2. Inner loop performs preprocessing, feature selection, class handling, and hyperparameter tuning on training data only.
3. The selected pipeline is refit on the full outer training set.
4. The outer test fold is evaluated exactly once.
5. Store predictions, selected features, hyperparameters, thresholds, and fold IDs.

## Primary classifier

SMO-trained SVM or closest reproducible implementation available in R.

Primary: linear kernel.  
Sensitivity: radial/RBF kernel if needed.

## Feature selectors

Primary:

1. t-test / ANOVA filter
2. Information Gain or Gain Ratio
3. SVM-RFE

Optional supplement:

- mRMR
- CFS
- Elastic net feature selection

## Class imbalance

Primary: class-weighted SVM.  
Sensitivity: SMOTE inside training folds only.  
Negative-control sensitivity: no imbalance correction.

## Metrics

- AUROC
- PR-AUC
- Balanced accuracy
- Matthews correlation coefficient
- Sensitivity
- Specificity
- Calibration/Brier score if probability outputs are reliable

## Leakage inflation metrics

```text
Delta_AUROC = AUROC_leaky - AUROC_guarded
Delta_PRAUC = PRAUC_leaky - PRAUC_guarded
Delta_MCC = MCC_leaky - MCC_guarded
```

## Feature stability

- Nogueira stability index
- Jaccard overlap
- Selection frequency
- Consensus feature list

## External validation rules

External validation cohorts must not be used for:

- feature selection
- threshold tuning
- hyperparameter tuning
- class balancing decisions
- batch correction that uses external labels or external cohort statistics in a way unavailable prospectively

## Reproducibility outputs

- code repository
- environment file
- random seeds
- fold assignments
- dataset accession table
- feature-selection frequency tables
- final evidence-audit table
- Zenodo archived release before submission

## Change-control rule

Any change to endpoint definition, dataset roles, feature selectors, or primary metrics after pilot analysis must be logged in this protocol.
