# v10.2 Final Consistency QA — Changelog

**Source:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_1_QA_patched.docx`
**Output:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_2_final_QA.docx`
**Date:** 2026-06-24
**Scope:** Text-only consistency patch. No analysis re-run; no result values changed; no references, code, or elastic-net results added; no BioTrust/founder/investor language introduced.

## Edits

1. **Results — "Cohort construction and endpoint harmonization" paragraph.**
   Replaced the prior wording that held GSE41998 out and stated the candidate cohorts were "not included in the results reported below." New wording states the cohort roles explicitly: GSE25055 = discovery/internal validation; GSE25065 = same-platform, same-study-family validation; GSE41998 (GPL571) = predeclared cross-platform transportability sensitivity cohort; GSE20194/GSE20271 = documented for audit transparency but held out owing to MDACC-lineage patient-overlap risk.

2. **Methods — "Study design overview."**
   Expanded the dataset sentence from GSE25055 + GSE25065 only to:
   "The present analysis uses GSE25055 for discovery/internal validation, GSE25065 for same-platform same-study-family validation, and GSE41998 as a predeclared cross-platform transportability sensitivity cohort."

3. **Scope statement.**
   Replaced "one discovery cohort and one same-study-family external cohort" with
   "one discovery cohort, one same-study-family validation cohort, and one cross-platform transportability sensitivity cohort."

4. **Results heading.**
   "External validation on GSE25065" → "Same-platform, same-study-family validation on GSE25065."

5. **Table 12 (Reproducible Omics Evidence Audit) — GSE25065 external-validation interpretation cell.**
   Replaced "Same-study-family validation shows a transportability drop; broader external transportability remains untested." with
   "Same-study-family validation shows a transportability drop; GSE41998 provides a separate cross-platform sensitivity check."
   The GSE41998 audit row is unchanged.

## Consistency sweep (verified)

- No "GSE41998 not used" remains.
- No "These cohorts are therefore not included in the results … below" remains for GSE41998.
- No "pre-registered" remains (0 occurrences).
- GSE25065 remains framed as same-platform, same-study-family validation.
- GSE41998 remains framed as a predeclared cross-platform transportability sensitivity cohort — not proof of leakage. Leakage and transportability remain distinct.
- The Limitations "remains untested" sentence (selector families / model classes) is a separate, accurate claim and was left unchanged.
- All GSE41998 numbers unchanged: N = 253; 69 pCR; 184 RD; AUROC 0.6638; PR-AUC 0.4353; balanced accuracy 0.5697; MCC 0.1845; sensitivity 0.2319; specificity 0.9076; confusion 16 TP / 53 FN / 167 TN / 17 FP; z-score sensitivity AUROC 0.6779 / PR-AUC 0.4300.
- All prior numbers unchanged: 0.7705, 0.7265, +0.044, 0.6078, 0.5409, 27 of 30, K-sweep values.
- No BioTrust/founder/investor language; no elastic-net results invented; no references added; no code edited.
- Table and figure numbering remain sequential; DOCX validates and renders (245 paragraphs, 13 tables).
