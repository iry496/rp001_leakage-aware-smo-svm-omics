# Reproducible Omics Evidence Audit - Notes (v1)

## Purpose
First integrated evidence audit for the leakage-aware SMO/SVM omics study.
It consolidates four completed stages - dataset audit, GSE25055 leaky-vs-guarded
pilot, feature-stability analysis, and GSE25065 external validation - into one
auditable table. No models were run here; all numbers are read from committed
outputs and recombined.

## Cohorts
- GSE25055 (discovery): RD=249, pCR=57, 4 NA excluded (N=310); pCR prevalence 18.6%.
- GSE25065 (external): RD=140, pCR=42, 16 NA excluded (N=198); pCR prevalence 23.1%.
- Same platform (GPL96), non-overlapping by design (GSE25066 split).

## Workflow sensitivity (naive vs guarded nested, GSE25055)
- Naive AUROC 0.7830 vs guarded nested AUROC 0.7032 -> whole-workflow contrast +0.0799.
- Naive PR-AUC 0.4257 vs guarded nested PR-AUC 0.3476 -> whole-workflow contrast +0.0781.
- The workflows also differ in tuning and resampling, so the contrasts are
  descriptive and not isolated causal estimates of feature-selection leakage.
  The guarded nested workflow improves the
  imbalance-aware metrics (balanced accuracy 0.5724, MCC 0.2138) over the leaky
  baseline (balanced accuracy 0.5525, MCC 0.1945).

## Feature stability (5 outer folds, K=100)
- 229 unique features selected across folds; 26 in all 5 folds (stable core);
  105 selected in exactly one fold (unstable tail).
- Mean Jaccard 0.3487, median Jaccard 0.3246, Nogueira stability 0.5128.
- Interpretation: stability is MODERATE - a small reproducible core coexists
  with a large unstable tail.

## External validation (GSE25065) and internal->external drop
- External AUROC 0.5633, PR-AUC 0.2986, balanced accuracy 0.5310, MCC 0.1013, sensitivity 0.1190, specificity 0.9429.
- Drop from guarded nested CV to external: AUROC +0.1399, PR-AUC +0.0491, balanced accuracy +0.0414,
  MCC +0.1125, sensitivity +0.0739, specificity +0.0090.
- Interpretation: external validation shows clear TRANSPORTABILITY LIMITS;
  discrimination and minority-class recall fall on the independent cohort,
  while specificity is roughly maintained.

## Cautious interpretation
- This is a METHODOLOGY / AUDIT result, NOT a clinical biomarker discovery.
- No clinical claims are made; pCR sensitivity is low internally and externally.
- The leaky pipeline inflates AUROC/PR-AUC; guarded nested validation gives an
  honest estimate and improves imbalance-aware metrics.
- Feature stability is moderate (stable core + unstable tail).
- External validation demonstrates real transportability limits.

## Limitations / unresolved risks
- Single same-platform external cohort (GSE25065).
- GSE41998 (cross-platform) and GSE20194/GSE20271 (overlap/de-dup) deliberately
  excluded pending harmonization and sample-level de-duplication.
- Feature-tail instability: 105 of 229 features selected only once.
- Low pCR sensitivity limits practical utility.

