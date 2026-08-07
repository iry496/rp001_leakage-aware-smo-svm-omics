# Reproducible Omics Evidence Audit

This repository supports the manuscript:

**A Reproducible Evidence-Audit Framework for Leakage, Feature Stability, and Transportability in Translational Omics Classification**

- Target journal: **Journal of Biomedical Informatics**
- Article type: **Research Paper**

This repository contains the submission-aligned analysis code, random seeds, fold assignments, selected-feature lists, software-environment files, generated figures and tables, supplementary materials, and preserved analysis outputs for manuscript version 1.2.2.

## Archived release

- Current submission-aligned release (v1.2.2): the version-specific DOI will be assigned when this release is archived.
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
                  figures, and evidence-audit dashboard
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

See [`CITATION.cff`](CITATION.cff). Until the v1.2.2 version DOI is assigned, cite the project using the concept DOI [10.5281/zenodo.21134085](https://doi.org/10.5281/zenodo.21134085).
