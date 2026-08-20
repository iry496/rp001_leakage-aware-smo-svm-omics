# Reproducible Omics Evidence Audit

This repository supports the manuscript:

**A Reproducible Evidence-Audit Framework for Leakage, Feature Stability, and Transportability in Translational Omics Classification**

- Target journal: **Journal of Biomedical Informatics**
- Article type: **Research Paper**

This repository contains the submission-aligned analysis code, random seeds, fold assignments, selected-feature lists, software-environment files, generated figures and tables, supplementary materials, and preserved analysis outputs for manuscript version 1.2.3.

## Archived release

- Current submission-aligned version DOI (v1.2.3): [10.5281/zenodo.21842032](https://doi.org/10.5281/zenodo.21842032)
- Previous five-author synchronization release DOI (v1.2.2): [10.5281/zenodo.21840696](https://doi.org/10.5281/zenodo.21840696)
- Previous synchronization release DOI (v1.2.1): [10.5281/zenodo.21840371](https://doi.org/10.5281/zenodo.21840371)
- Previous submission-aligned version DOI (v1.2.0): [10.5281/zenodo.21834590](https://doi.org/10.5281/zenodo.21834590)
- Prior version DOI (v1.1.0): [10.5281/zenodo.21134086](https://doi.org/10.5281/zenodo.21134086)
- Concept DOI (all versions): [10.5281/zenodo.21134085](https://doi.org/10.5281/zenodo.21134085)

## Summary

High-dimensional omics classifiers can appear credible when data leakage, unstable feature selection, class-imbalance behavior, and weak external transportability remain hidden. This work presents a reproducible **evidence-audit framework** — the Reproducible Omics Evidence Audit — that integrates leakage sensitivity, a label-permutation negative control, guarded (nested) validation, feature-selection stability, external and cross-platform transportability, class-imbalance behavior, reproducibility artifacts, and explicit red-flag triggers into a single reusable reporting instrument. Public breast-cancer neoadjuvant-chemotherapy cohorts (pCR vs. residual disease) serve as a high-dimensional stress test; an established linear SVM is used as a transparent workhorse, not as a methodological advance.

## Submission-aligned results

| Audit component | Result |
| --- | --- |
| Naive workflow, GSE25055 | AUROC 0.7830; PR-AUC 0.4257 |
| Guarded nested workflow, GSE25055 | AUROC 0.7032; PR-AUC 0.3476 |
| Descriptive whole-workflow contrast | ΔAUROC 0.0799 (95% CI 0.0080–0.1514; p=0.026); ΔPR-AUC 0.0781 (−0.0423–0.1755; p=0.192) |
| Fixed-orientation permutation control | Naive null AUROC mean 0.8778; guarded null mean 0.4917 |
| Thirty-seed robustness | AUROC contrast positive in 28/30 seeds; PR-AUC contrast positive in 22/30 |
| Feature stability | Nogueira 0.5128; mean Jaccard 0.3487; 26 probes in all five anchor folds |
| GSE25065 | AUROC 0.5633; PR-AUC 0.2986 |
| GSE41998 primary | AUROC 0.6922; PR-AUC 0.4634 (253 explicit Yes/No cases) |

The naive and guarded arms differ in feature-selection placement, cost policy, and resampling. Their difference is therefore a descriptive whole-workflow contrast, not an isolated causal estimate of leakage.

## Current working analysis extensions

These post-v1.2.3 extensions are reproducible working-branch analyses and are
not part of the archived submission-aligned release above.

| Extension | Result |
| --- | --- |
| Thirty-repeat matched ablation | Leaky mean AUROC 0.7870; guarded 0.7256; paired mean difference 0.0615; positive in 30/30 repeats; empirical 2.5th–97.5th percentile split interval 0.0084–0.1201 |
| Matched PR-AUC | Leaky mean 0.4237; guarded 0.3733; paired mean difference 0.0505; positive in 24/30 repeats |
| Unified GSE25065 projection | AUROC 0.5633 (conditional 95% bootstrap interval 0.4643–0.6684); calibration slope 0.1480; sensitivity/specificity 0.1190/0.9429 |
| Unified GSE41998 primary projection | AUROC 0.6922 (0.6146–0.7599); calibration slope 0.6166; sensitivity/specificity 0.7971/0.5000 |

The matched ablation uses identical outer and inner folds, tuning policy, cost
grid, scaling, class weighting, model family, probability procedure, final-fit
random-number seed, 0.5 threshold, and estimators in both arms. Only supervised
feature-selection placement changes. Its split-distribution interval describes
robustness across correlated repeated partitions and is not an
independence-based confidence interval.

The unified external analysis fits one GSE25055 model once and reuses that exact
model, discovery-derived scaler, and 0.5 threshold for both primary external
projections. Its bootstrap intervals are conditional on the frozen model and
stored external predictions; they do not include discovery-workflow refitting
uncertainty.

## Cohorts

| Cohort | Platform | Role |
| --- | --- | --- |
| GSE25055 | Affymetrix HG-U133A (GPL96) | Discovery / internal validation |
| GSE25065 | Affymetrix HG-U133A (GPL96) | Same-platform, same-study-family validation |
| GSE41998 | Affymetrix HG-U133A 2.0 (GPL571) | Cross-platform transportability sensitivity |
| GSE20194 / GSE20271 | — | Documented but held out (MDACC-lineage patient-overlap risk) |

For GSE41998, the live Series Matrix contains 69 explicit `Yes`, 184 explicit `No`, 20 literal `0`, and 6 missing values in the `pcr` field. The zeros remain a distinct category in `pcrrcb1`; the analysis maps Yes→pCR and No→RD and excludes zeros/missing values without assigning an undocumented clinical meaning.

## Repository structure

```
R/                Reusable R functions (feature selection, preprocessing, model, metrics)
scripts/          Analysis scripts: dataset audit, leaky baseline, guarded nested pipeline,
                  feature stability, external validation (GSE25065, GSE41998), evidence-audit
                  table, bootstrap CIs, permutation control, repeated nested CV, selector K-sweep,
                  matched one-factor ablation, unified external validation, figures, and
                  evidence-audit dashboard
notebooks/        Quarto notebooks
data_accessions/  GEO accession registry
processed_data/   Derived matrices (large files not committed)
raw_data/         Local raw downloads (not committed)
results/          Analysis outputs: metrics, predictions, selected features, fold assignments,
                  permutation null distributions, bootstrap CIs
figures/          Generated figures
tables/           Generated tables (dataset audit, pipeline comparison, evidence audit)
environment/      Package/version environment files
manuscript/       Manuscript files
supplementary/    Supplementary files
```

## Leakage-control rules

| Step | Guarded implementation |
|---|---|
| Scaling | Fit on training fold only |
| Supervised feature selection | Training fold only |
| Hyperparameter tuning | Inner cross-validation loop only |
| Class weighting | Training fold only |
| External validation | Frozen preprocessing, feature set, model, and threshold |

## Reproduce the submission outputs

Use R 4.4.x with the packages declared in `environment/packages.R`. From the
repository root, run the numbered scripts in order. The submission-critical
refresh sequence is:

```sh
Rscript scripts/07_bootstrap_ci.R
Rscript scripts/10_build_evidence_audit_artifact.R
Rscript scripts/11_plot_evidence_audit_dashboard.R
```

Run the working analysis extensions with:

```sh
Rscript scripts/16_matched_ablation_gse25055.R full 30
Rscript scripts/17_external_validation_unified.R
```

`scripts/07_bootstrap_ci.R` is the sole canonical uncertainty implementation
because it uses the manuscript-declared pROC and PRROC estimators. The Python
script under `scripts/checks/` writes separate approximation files only. Raw
expression matrices are retrieved from GEO and are not redistributed; the
accessions, committed predictions, folds, selected features, expected tables,
and figures provide an auditable reproduction path.

## License

This repository is dual-licensed:

- **Code** (`R/`, `scripts/`, `notebooks/`) under the **MIT License** — see [`LICENSE`](LICENSE).
- **Manuscript, figures, tables, and derived data content** under **Creative Commons Attribution 4.0 International (CC-BY-4.0)** — see [`LICENSE-CONTENT`](LICENSE-CONTENT).

## How to cite

See [`CITATION.cff`](CITATION.cff). Cite the submission-aligned v1.2.3 archive using DOI [10.5281/zenodo.21842032](https://doi.org/10.5281/zenodo.21842032).
