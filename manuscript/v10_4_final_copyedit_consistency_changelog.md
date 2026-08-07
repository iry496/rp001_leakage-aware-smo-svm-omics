# v10.4 — Final Desk-Review Copyedit & Consistency — Changelog

**Source:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_3_desk_review_blocker_cleanup.docx`
**Output:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_4_final_copyedit_consistency.docx`
**Date:** 2026-06-24
**Scope:** Copyedit/consistency only. No analysis; no new models; no elastic-net results; no new references; no code edits; no result values changed; no BioTrust/founder/investor language.

## Change made (exactly one)

### Task 1 — Capitalization / grammar
- Sentence-start fix: "This study does not overstate algorithmic novelty. **linear SVM is not new**; the novelty is the audit architecture around it." → "… **The linear SVM is not new**; …" (paragraph beginning with the "Models and methods / novelty" discussion).

That is the only textual change in the document. A whole-document diff confirms exactly one paragraph differs from v10.3, and all numeric tokens are byte-identical between v10.3 and v10.4.

### Other Task 1 items checked
- "linear SVM is not treated…" already reads "The linear SVM is not treated here as a new algorithm." (fixed in v10.3) — unchanged.
- Full scan for other sentence-start lowercase model/pipeline phrases ("linear SVM", "naive", "guarded", "leaky", "SMO", etc.): none found other than the para 44 case above. (Intentional lowercase starts left as-is: the "pCR =" table legend, the "nestedcv:" package-name article title, and DOI lines in the reference list.)

## Items reviewed and intentionally left unchanged

### Task 2 — Placeholder scan
- "Affiliation to be finalized" / "Affiliation not provided": **none** (Author 3, Paul Tan, is listed as "Harvard University Extension School"; no affiliation placeholder remains).
- "TBD", "to be finalized", "placeholder": **none**.
- "planned" (implying completed work is still planned): **none**. The only forward-looking "to be reported once executed" text concerns the elastic-net comparator, which is genuinely future work (correct).
- **Remaining submission-blocking admin items (not auto-completable — require author input; nothing invented):**
  - "Competing interests will be declared by all authors prior to submission."
  - "Funding sources will be declared prior to submission."
  - "Author contributions will be finalized prior to submission."
  - "Acknowledgements will be added prior to submission."
  These four standard pre-submission declarations remain **submission-blocking items** for the authors to complete.

### Task 3 — SMO/SVM wording
- General prose already uses "linear SVM." Solver detail ("SMO-trained linear SVM" / "trained by sequential minimal optimization") appears only where appropriate (Abstract method intro, Scope statement, Methods), and the Platt (1999) SMO reference is retained. The manuscript does not read like an SMO/SVM algorithm paper. No change needed.

### Task 4 — Competitor positioning
- nestedcv and OmicSelector are described as existing tools providing leakage-resistant machinery ("Remaining methodological gap" paragraph and the contribution paragraph: "This work does not replace nestedcv, stabm, OmicSelector … the present contribution is the empirical audit such machinery enables"). The manuscript is framed as a reproducible empirical stress-test audit and evidence-reporting artifact, not a competing package. No unverified claims about OmicSelector. No change needed.

## Task 5 — Internal reference check (verified)
- Figure references sequential and complete: Figures 1–8.
- Table references sequential and complete: Tables 1–12. Table 12 / Reproducible Omics Evidence Audit references correct.
- GSE25065 remains "same-platform, same-study-family validation."
- GSE41998 remains "cross-platform transportability sensitivity."
- No "GSE41998 not used" and no "not included below / not included in the results" language remains.

## Task 6 — Scientific integrity (verified)
- No numerical result changed (numeric tokens identical to v10.3); all key values retained: 0.7705, 0.4016, 0.7265, 0.3653, ΔAUROC +0.044 (−0.013 to +0.103), 27 of 30, +0.054, 0.5409, 0.6078, 0.3060, 0.6638, 0.4353, GSE25065 CI 0.506–0.700, N = 253.
- Elastic net remains future work only (paras "candidate supplementary analyses" and the Future-work elastic-net comparator).
- No biomarker discovery claim; recurrent probes described as stability-ranked candidates.
- No patient-level overclaiming from the 30-seed Wilcoxon result (framed as seed/fold-level robustness, not patient-level inference).
- Leakage and transportability remain distinct.
- No "pre-registered"; no BioTrust/founder/investor language.
- DOCX validates and renders (245 paragraphs, 13 tables).
