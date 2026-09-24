# v10.3 — Desk-Review Blocker Cleanup — Changelog

**Source:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_2_final_QA.docx`
**Output:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_3_desk_review_blocker_cleanup.docx`
**Date:** 2026-06-24
**Scope:** Mechanical desk-review blocker cleanup only. No analysis re-run; no result values changed (values were de-duplicated or relocated, never altered); no new references; no code edits; no elastic-net results; no BioTrust/founder/investor language; scientific argument unchanged.

## Changes made

### Task 3 — Linear SVM typo (fixed)
- Introduction paragraph that began "linear SVM is not treated here as a new algorithm." now reads "**The** linear SVM is not treated here as a new algorithm."

### Task 5 — Abstract density (lightly shortened: 2,511 → 2,407 characters)
Removed duplicated / non-essential numeric detail from the Abstract; every value removed already appears in Results and/or Tables (relocation, not deletion). No value was changed.
- Removed the per-estimate AUROC bootstrap CIs from the Abstract — "AUROC 0.7705 (95% bootstrap CI 0.711–0.824)" → "AUROC 0.7705"; "AUROC 0.7265 (0.653–0.791)" → "AUROC 0.7265". (Both CIs remain in Table 4 and in the Results CI paragraph.)
- Removed the duplicated paired-bootstrap p-value from the Abstract leakage-gap parenthetical — "(ΔAUROC +0.044, 95% CI −0.013 to +0.103; paired bootstrap p = 0.13)" → "(ΔAUROC +0.044, 95% CI −0.013 to +0.103)". The "CI includes zero" detail is retained. (p = 0.13 remains in Results.)
- Relocated the GSE25065 external-validation AUROC CI from the Abstract to Results: removed "(0.506–0.700)" from the Abstract and added "(95% CI 0.506–0.700)" to the GSE25065 Results sentence ("the model achieved AUROC 0.6078 (95% CI 0.506–0.700), PR-AUC 0.3060 …"). The value is unchanged; it now lives where the other validation metrics are reported.
- Trimmed the duplicated platform/sample parenthetical from the Abstract — "the independent GSE41998 cohort (GPL571; N = 253) yielded AUROC 0.6638" → "the independent GSE41998 cohort yielded AUROC 0.6638". (GPL571 and N = 253 remain in Methods, Results, and Table 12.)

**Essential numbers retained in the Abstract:** GSE25055 sample count and class balance (306; 57 pCR, 249 RD); leaky vs guarded AUROC/PR-AUC (0.7705 / 0.4016 vs 0.7265 / 0.3653); ΔAUROC +0.044 with 95% CI −0.013 to +0.103 (CI includes zero); 1000-permutation negative control; 30-seed directionality (27 of 30; median +0.054); feature stability (Nogueira 0.5409); GSE25065 AUROC 0.6078; GSE41998 AUROC 0.6638; contribution statement.

## Items reviewed and intentionally left unchanged

### Task 1 — Author / admin fields
- **Affiliations:** No affiliation placeholder remains in the body. Author 3 (Paul Tan) now shows "Harvard University Extension School"; no "Affiliation to be finalized / not provided" text exists. Nothing invented.
- **Remaining submission-blocking admin items (NOT auto-completable — flagged here for the authors):** the back-matter still contains four standard pre-submission declarations that require author input and were left as-is because their content cannot be truthfully supplied by an editor:
  - "Competing interests will be declared by all authors prior to submission."
  - "Funding sources will be declared prior to submission."
  - "Author contributions will be finalized prior to submission."
  - "Acknowledgements will be added prior to submission."
  These remain **submission-blocking admin items** to be filled in by the authors before submission.

### Task 2 — Truncated / incomplete sentences
- Full scan for ellipses ("…", "..."), abrupt endings, and dangling citations across the flagged sections (Remaining methodological gap, Direct Competitor Positioning, nestedcv, OmicSelector, GSE41998, Evidence Audit) and the whole body. **No truncated or incomplete sentences found.** Nothing changed.

### Task 4 — SMO/SVM wording
- General prose already uses "linear SVM." Solver detail ("SMO-trained linear SVM" / "trained by sequential minimal optimization") appears only where appropriate: the Abstract method intro, the Scope statement (explicitly stating SMO is only the solver, not a new algorithm), and Methods. The Platt (1999) SMO reference is retained. No change needed; nothing reads like an SMO/SVM algorithm paper.

### Task 6 — Competitor positioning
- nestedcv and OmicSelector are already framed as existing tools providing reusable leakage-resistant machinery (Remaining-methodological-gap paragraph and the contribution paragraph: "This work does not replace nestedcv, stabm, OmicSelector … the present contribution is the empirical audit such machinery enables"). The manuscript is already framed as an empirical stress-test audit and standardized evidence-reporting artifact, not a competing package. No unverified claims about OmicSelector were added. No change needed.

## Task 7 — Consistency checks (verified)
- No "pre-registered" (0). No "GSE41998 not used" (0). No BioTrust/founder/investor language (0 each).
- No GSE41998 future/planned wording; future-work text refers only to further cohorts beyond GSE41998.
- GSE25065 remains "same-platform, same-study-family validation"; GSE41998 remains "cross-platform transportability sensitivity"; leakage and transportability remain distinct.
- Figures 1–8 and Tables 1–12 remain sequential.
- All key result values unchanged: 0.7705, 0.4016, 0.7265, 0.3653, +0.044 (−0.013 to +0.103), 27 of 30, +0.054, 0.5409, 0.6078, 0.3060, 0.6638, 0.4353, N = 253, plus the relocated GSE25065 CI 0.506–0.700.
- DOCX validates and renders (245 paragraphs, 13 tables).
