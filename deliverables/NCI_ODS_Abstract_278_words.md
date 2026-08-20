# NCI ODS 2026 Poster Abstract

## Category

Poster Abstract

## Proposed title

From Shared Cancer Data to Auditable Evidence: A FAIR Audit of Breast-Cancer pCR Prediction

## Keywords

Cancer data sharing; FAIR reuse; breast cancer transcriptomics; information leakage; transportability

## Abstract Summary — 278 words

Publicly shared cancer transcriptomic data can accelerate discovery through independent reuse and cross-cohort validation, but availability alone does not make resulting evidence reliable. We developed the Reproducible Omics Evidence Audit (ROEA), a FAIR-oriented workflow linking accession-level provenance and label harmonization to seeds, folds, selected features, predictions, software versions, and human- and machine-readable outputs. We applied ROEA to three Gene Expression Omnibus breast-cancer cohorts with pathological complete response after neoadjuvant therapy as the endpoint: GSE25055 (n=306), GSE25065 (n=182), and GSE41998 (n=253).

A 30-repeat matched ablation held folds, tuning, scaling, class weighting, model family, threshold, and estimators constant while changing only supervised feature-selection placement. Selecting features before cross-validation increased mean AUROC from 0.726 to 0.787 (difference 0.061; empirical split interval 0.008–0.120; positive in 30/30 repeats); the PR-AUC difference was smaller and less consistent (0.051; positive in 24/30 repeats). This ablation was independently reproduced across R versions from the supplied cached matrix. In separate negative controls using earlier unmatched complete workflows, the naive permutation-null AUROC averaged 0.878—higher than its own real-label AUROC of 0.783—whereas the guarded null centered near chance (0.492). Guarded selection stability was moderate (Nogueira index 0.51); recurrent probes are audit outputs, not validated biomarkers.

One GSE25055 model was fit once and projected without outcome-informed retuning. External AUROCs were 0.563 in GSE25065 (conditional 95% bootstrap interval, 0.464–0.668) and 0.692 in GSE41998 (0.615–0.760). The same 0.50 threshold yielded sensitivity/specificity of 0.119/0.943 and 0.797/0.500. Both projections had negative Brier skill and calibration slopes below 1, showing that discrimination did not ensure transportable operating behavior.

Public sharing enabled independent reuse, workflow-specific leakage testing, and frozen cross-cohort evaluation. ROEA provides a reusable, versioned audit pattern for provenance-controlled, leakage-tested, stability-audited, and externally challenged evidence.

No DOI is included in this submission version.
