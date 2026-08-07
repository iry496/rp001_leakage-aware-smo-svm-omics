# Changelog

## 1.2.3 — 2026-08-07

Final minor synchronization release for the Journal of Biomedical Informatics submission.

- Replaced the stale Section 3.11 statement of uniformly low pCR sensitivity with cohort- and scaling-dependent operating-point behavior for GSE25065 and GSE41998.
- Aligned the Declaration of Interests with Iris Yang's co-corresponding-author role.
- Updated the retained protocol provenance note to identify v1.2.2 as the superseding release.
- Recorded completion of all five final-author approvals and the written consent for the authorship change in the cover letter.
- Regenerated the manuscript, editorial documents, QA proof, and submission package without changing scientific results.
- Archived the final minor synchronization release at Zenodo DOI 10.5281/zenodo.21842032.

## 1.2.2 — 2026-08-07

Bounded authorship-synchronization release for the Journal of Biomedical Informatics submission.

- Standardized the five-author order as Iris Yang, Paul Tan, Jung Chen, Jewel Wang, and Chung-I Huang.
- Confirmed Jung Chen as the legal publication name and University of Chicago as the affiliation.
- Retained Iris Yang as first and corresponding author and Chung-I Huang as co-corresponding author.
- Assigned Jung Chen the same CRediT roles as Paul Tan: Validation; Writing – Review & Editing.
- Removed the former sixth author and associated affiliation from the current manuscript, editorial documents, repository metadata, and release package.
- Regenerated the manuscript, editorial documents, QA proof, and submission package without changing scientific results.
- Archived the five-author synchronization release at Zenodo DOI 10.5281/zenodo.21840696.

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
