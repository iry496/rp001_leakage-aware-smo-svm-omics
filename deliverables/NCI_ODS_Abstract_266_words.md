# NCI ODS 2026 Poster Abstract

## Category

Poster Abstract

## Proposed title

From Shared Cancer Data to Auditable Evidence: A FAIR Audit of Breast-Cancer pCR Prediction

## Keywords

Cancer data sharing; FAIR reuse; breast cancer transcriptomics; information leakage; transportability

## Abstract Summary — 266 words

Publicly shared cancer transcriptomic data can accelerate discovery through independent reuse and cross-cohort validation, but availability alone does not make resulting evidence reliable. We developed the Reproducible Omics Evidence Audit (ROEA), a FAIR-oriented workflow linking accession-level provenance and label-harmonization decisions to seeds, folds, selected features, predictions, software versions, and human- and machine-readable outputs. We applied ROEA to three Gene Expression Omnibus breast-cancer cohorts with pathological complete response after neoadjuvant therapy as the endpoint: GSE25055 (n=306), GSE25065 (n=182), and GSE41998 (n=253).

A 30-repeat matched ablation held outer and inner folds, tuning, scaling, class weighting, model family, threshold, and estimators constant while changing only supervised feature-selection placement. Selecting features before cross-validation increased mean AUROC from 0.726 to 0.787 (paired mean difference, 0.061; positive in 30/30 repeats; empirical 2.5th–97.5th percentile split interval, 0.008–0.120). In separate complete-workflow negative controls with 1,000 label permutations, the naive workflow’s null AUROC averaged 0.878, whereas the guarded null centered near chance (0.492). Guarded selection stability was moderate (Nogueira index, 0.51); recurrent probes are audit outputs, not validated biomarkers.

One GSE25055 model was fit once and projected without outcome-informed retuning. External AUROCs were 0.563 in GSE25065 (conditional 95% bootstrap interval, 0.464–0.668) and 0.692 in GSE41998 (0.615–0.760). The identical 0.50 threshold yielded sensitivity/specificity of 0.119/0.943 and 0.797/0.500, respectively. Both discovery-scaled projections had negative Brier skill, and calibration slopes were 0.148 and 0.617, demonstrating that discrimination did not ensure transportable operating behavior.

Public sharing enabled independent secondary analysis, workflow-specific leakage testing, matched ablation, and frozen cross-cohort evaluation. ROEA provides a reusable, versioned audit pattern for turning shared cancer data into provenance-controlled, leakage-tested, stability-audited, and externally challenged evidence.

No DOI is included in this submission version.
