# Phase 1 Completion Summary

_Status snapshot for the BioData Mining methodology manuscript_
**"An Audit of Data Leakage and Feature Stability in Transcriptomic Biomarker Classification Using Nested SMO/SVM and External Validation."**

This summary reflects only material already merged into `main` (PRs #19–#22). It introduces no new analyses, manuscript edits, or code, and changes no result values.

---

## 1. What Phase 1 completed

Phase 1 ("Foundations") moved the project from a draft with pilot results to a **Prof. Wong / Paul review-ready package**: statistically quantified results, a clean completed-study manuscript, and corrected novelty framing. Concretely, Phase 1 delivered:

- **Statistical uncertainty** for every headline metric (bootstrap confidence intervals) plus a significance test for the leakage gap (paired bootstrap + DeLong).
- **Manuscript hygiene** — converted the draft from proposal/planning language to a completed-study report.
- **Novelty reframe** — positioned the contribution as an integrated empirical audit rather than a new algorithm.

The revision roadmap and the GitHub issue/label system that organize Phases 1–4 were also established (issues #7–#18).

## 2. PRs merged (#19–#22)

| PR | Branch | Merge commit | Summary |
|----|--------|--------------|---------|
| #19 | `docs/revision-roadmap` | `b829a9a` | Added `docs/revision_roadmap_after_ai_review.md` (the phased plan). |
| #20 | `analysis/bootstrap-ci` | `70dbd95` | Added bootstrap CIs, paired-bootstrap and DeLong tests, R script + Python cross-check + notebook. |
| #21 | `manuscript/hygiene-tense-cleanup` | `e978432` | Cleaned tense/notes/placeholders; fixed an orphan reference. |
| #22 | `manuscript/novelty-reframe` | `cc58664` | Reframed novelty as an integrated audit; positioned vs existing tools. |

All four branches remain on the remote (not deleted).

## 3. Current manuscript version and title

- **Latest version:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v2_phase1_novelty_reframe.docx` (built on the hygiene-cleaned `…v2_phase1_clean.docx`).
- **Title:** *An Audit of Data Leakage and Feature Stability in Transcriptomic Biomarker Classification Using Nested SMO/SVM and External Validation.*
- **Framing:** methodology / audit paper; explicitly **not** a clinical biomarker discovery study.

## 4. Key results from the leakage audit (GSE25055)

Leaky baseline vs. guarded nested SMO/SVM on the discovery cohort (306 usable samples; 57 pCR, 249 RD):

| Metric | Leaky | Guarded nested | Gap (leaky − guarded) |
|--------|------:|---------------:|----------------------:|
| AUROC | 0.7705 | 0.7265 | +0.0440 |
| PR-AUC | 0.4020 | 0.3656 | +0.0363 |
| Balanced accuracy | 0.5546 | 0.5792 | −0.0246 |
| MCC | 0.2063 | 0.2250 | −0.0187 |
| Sensitivity (pCR) | 0.1333 | 0.2105 | −0.0772 |
| Specificity (RD) | 0.9759 | 0.9478 | +0.0281 |

The leaky pipeline inflated the threshold-independent discrimination metrics (AUROC, PR-AUC). Operating-point metrics moved in the opposite direction but are threshold-dependent under class imbalance and are interpreted cautiously; minority-class sensitivity remained low under both pipelines.

## 5. Key feature-stability results

Across the five outer folds of the guarded nested pipeline (top 100 probes per fold from a 22,283-probe universe):

- **222** unique probes selected at least once.
- **28** probes selected in all five folds (stable core); **102** selected in only one fold (unstable tail).
- Mean pairwise Jaccard **0.3734**, median **0.3699**; Nogueira stability index **0.5409**.
- Interpretation: **moderate** stability — a small reproducible core coexists with a large unstable tail. Selected features are audit outputs, not validated markers.

## 6. Key external-validation results (GSE25065)

Frozen guarded model (linear SVM, cost 0.25, 100 discovery-selected probes, all aligned) projected once onto GSE25065 (182 usable; 42 pCR, 140 RD):

- AUROC **0.6078**, PR-AUC **0.3060**, balanced accuracy **0.5381**, MCC **0.1347**, sensitivity **0.1190**, specificity **0.9571** (confusion: 5 TP / 6 FP / 134 TN / 37 FN).
- **Transportability drop** vs guarded nested CV: AUROC −0.1187, PR-AUC −0.0596, balanced accuracy −0.0411, MCC −0.0903, sensitivity −0.0915; specificity essentially unchanged (+0.0094).
- This decline reflects cross-cohort transportability limits, distinct from the internal leakage gap.

## 7. What bootstrap CI added (#20)

- **95% confidence intervals** (stratified percentile bootstrap, B = 2000, seed 20260620) for AUROC, PR-AUC, balanced accuracy, MCC, sensitivity, and specificity — for the leaky, guarded, and external results.
- **Significance test for the leakage gap**: paired stratified bootstrap and DeLong.
- **Key finding (feeds claim calibration):** the leakage gap is **positive but not statistically significant** in this single cohort.
  - ΔAUROC +0.0440, 95% CI **[−0.0125, +0.1028]**, bootstrap p ≈ 0.13.
  - ΔPR-AUC +0.0362, 95% CI **[−0.0483, +0.1247]**, bootstrap p ≈ 0.38.
  - DeLong per-repeat p-values 0.07–0.51 (median 0.18).
  - All CIs are wide, driven by the small pCR minority (57/306 internally, 42/182 externally).
- **Files:** `scripts/07_bootstrap_ci.R` (canonical, pROC + PRROC), `scripts/checks/07_bootstrap_ci_crosscheck.py`, `notebooks/07_bootstrap_ci.qmd`, `tables/uncertainty/{bootstrap_ci.csv, delta_auroc_prauc_ci.csv}`, `results/uncertainty/{delong_tests.csv, bootstrap_ci_notes.md}`. No new modeling — resampling of stored out-of-fold predictions.

## 8. What hygiene cleanup changed (#21)

Text only — no result values changed (187 paragraphs, 9 tables preserved):

- Proposal/future-tense ("proposed manuscript", "will report", "planned", "must avoid") → completed-study tense.
- Removed self-instructions / internal editorial notes (e.g., "This manuscript should avoid…", the Table 2 "preserve in the cover letter / response to reviewers" note, the BioTrust placement instruction, the ethics "should be reviewed…" note).
- Cleaned "To be finalized" placeholders (Competing interests, Funding, Authors' contributions, Acknowledgements) and removed informal first-name planning.
- **Fixed an orphan reference:** Berrar & Flach (2012) was listed but never cited → now cited in the metrics paragraph; no reference deleted.

## 9. What novelty reframe changed (#22)

Text/positioning only — no result values changed; only the qualitative competitor table (Table 2) was edited:

- States explicitly that the SMO/SVM classifier, nested cross-validation, and feature-stability metrics are **established methods, not the novelty**.
- Frames the contribution as the **integrated empirical audit**: leakage gap + feature stability + external validation + evidence-audit reporting.
- Added explicit positioning vs **nestedcv** (Lewis et al., 2023), **OmicSelector** (named as an automated omics feature-selection / deep-learning environment), and the data-leakage literature (Ambroise & McLachlan 2002; Kapoor & Narayanan 2023); contrasts pCR/RD signature studies.
- Added one OmicSelector row to the Direct Competitor Positioning table.
- Reinforced that the work is not a clinical biomarker discovery study.

## 10. Remaining limitations before submission

1. The leakage gap is **not statistically significant** in this single cohort; substantiating the leakage claim needs the permutation control and repeated nested CV.
2. **Single same-platform external cohort** (GSE25065); cross-platform transportability is untested.
3. **Single feature selector** (univariate t-test top-K) and **single classifier** (linear SMO/SVM).
4. **Low pCR sensitivity** internally and externally; limited practical utility for identifying responders.
5. **Moderate feature stability** with a large unstable tail (102 of 222 probes selected once).
6. **Reproducibility release** (environment lockfile, seeds/fold assignments, Zenodo DOI) not yet prepared; Declarations (competing interests, funding, contributions, acknowledgements) still hold clean placeholders to finalize.
7. The committed bootstrap CIs were generated by the Python cross-check because R is not runnable in the working environment; the canonical `scripts/07_bootstrap_ci.R` should be run in RStudio to confirm reproducibility (point estimates and DeLong are deterministic; CIs are reproducible within each language).

## 11. Recommended Phase 2 tasks

Per the revision roadmap (`docs/revision_roadmap_after_ai_review.md`):

- **#11 `analysis/permutation-control`** (P1, new modeling) — label-shuffle negative control under both pipelines; report null distributions and empirical p-values.
- **#12 `analysis/repeated-nested-cv`** (P1, new modeling) — repeat nested CV across 20–30 seeds; report leakage-gap and feature-stability distributions.
- Then finalize **#9 `manuscript/claim-calibration`** using the CI, permutation, and repeated-CV evidence.

(Later phases: #13 `analysis/gse41998-external` and #14 `analysis/comparator-elasticnet`; then #15 `figure/evidence-audit-dashboard` and #16 `repo/reproducibility-release`; optional #17 `analysis/selector-and-k-sweep`, #18 `analysis/calibration-threshold`.)

## 12. Questions for Prof. Wong (modeling / validation rigor)

1. For the leakage-gap significance test, is per-repeat DeLong plus a paired stratified bootstrap the right approach for the repeated-CV leaky baseline, or do you prefer a different test?
2. How many seeds for repeated nested CV (#12) — 20 or 30 — and any preferred outer/inner fold counts?
3. Is the univariate t-test top-K acceptable as the **primary** selector for this version, with the selector/K sweep (#17) as a sensitivity analysis, or should additional selectors (Information Gain, SVM-RFE, mRMR) be in the main analysis?
4. Is linear SMO/SVM sufficient as the sole classifier for the main paper, or should an RBF and/or elastic-net comparator (#14) be promoted from supplementary to main?
5. For GSE41998 (#13), please confirm the leakage-safe cross-platform design (no joint normalization, no global ComBat; training-derived probe/gene mapping; frozen GSE25055 model) and the GSE41998 label-coding resolution.

## 13. Questions for Paul (evidence-audit framing / translational relevance)

1. Does the "audit-not-algorithm" framing and the Reproducible Omics Evidence Audit table land well for the intended audience?
2. Is the evidence-audit dashboard figure (#15) worth prioritizing for submission, and which panels (leakage gap with CIs, stability, external drop, reproducibility status) matter most?
3. Is there BioTrust / translational framing you want included — kept in Discussion / future work and out of the title and abstract — and where should the line be?
4. Please confirm competing-interest disclosures and any acknowledgements/funding to record; the Declarations section currently holds clean placeholders pending finalization.
