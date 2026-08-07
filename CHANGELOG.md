# Changelog

## 1.2.1 — 2026-08-07

Bounded synchronization release for the Journal of Biomedical Informatics submission.

- Restored the manuscript-declared R/PRROC uncertainty workflow as the sole canonical source and regenerated all downstream evidence-audit outputs.
- Isolated the approximate Python cross-check so it can no longer overwrite submission tables.
- Added selector-K sensitivity results to Supplementary File S2.
- Clarified operating-point aggregation, GSE41998 thresholding, supplementary cross-references, archive metadata, and the Hatzis reference.
- Regenerated the manuscript, QA proof, S5 dashboard, repository metadata, and Drive submission package.
- Archived the synchronized release at Zenodo DOI 10.5281/zenodo.21840371.

## 1.2.0 — 2026-08-07

Submission-matched release for the Journal of Biomedical Informatics manuscript.

Zenodo version DOI: https://doi.org/10.5281/zenodo.21834590

- Fixed pROC direction a priori so resampling and permutation analyses cannot select ROC orientation from observed labels.
- Sorted GEO samples by accession before seeded fold creation, making the release invariant to download order.
- Prevented integer overflow in Matthews correlation coefficient counts.
- Regenerated the naive and guarded predictions, metrics, feature lists, uncertainty intervals, 1000-permutation control, 30-seed robustness results, selector-K sensitivity, figures, and evidence-audit artifacts.
- Reframed the naive-versus-guarded result as a descriptive whole-workflow contrast because the arms also differ in tuning and resampling.
- Re-executed GSE25065 and GSE41998 frozen-model validation; added GSE41998 bootstrap confidence intervals.
- Verified the GSE41998 binary endpoint from the live Series Matrix: Yes→pCR (69), No→RD (184), with 20 literal zeros and 6 missing values excluded.
- Synchronized six-author order and affiliations; Jung Chen has the same CRediT roles as Paul Tan (Validation; Writing – Review & Editing).
- Removed superseded manuscript drafts, internal review packets, and pre-correction result files from the release. They remain available in Git history.
- Updated the graphical abstract and supplementary materials to the same submission-aligned results.
