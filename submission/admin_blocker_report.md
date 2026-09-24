# Admin Blocker Report

**Scope:** Administrative / manuscript-support cleanup only. No model analysis was run; no scientific result values were changed; no classifiers, feature selectors, or elastic-net results were added; no author declarations, affiliations, funding, or contributions were invented.

**Source manuscript:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_5_table_figure_readability.docx`
**Updated manuscript (admin only):** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_6_admin_blocker_cleanup.docx`
**Date:** 2026-06-24

## Legend
- **Claude can fill:** safe to complete without inventing author facts — done.
- **Author input required:** cannot be completed without the authors — flagged, not invented.

---

## 1. Competing interests
- **Current text (manuscript v10.5 body):** "Competing interests will be declared by all authors prior to submission."
- **Why it blocks desk review:** an explicit "will be declared" placeholder signals an unfinished draft; BMC requires a real competing-interests statement at submission.
- **Author confirmation needed:** each of the five authors must declare any financial/non-financial competing interests, or confirm none.
- **Proposed final text (if all confirm none):** "The authors declare that they have no competing interests." (Use only if every author confirms.)
- **Who fills:** **Author input required.** Placeholder **removed from manuscript body** (v10.6); tracked in `competing_interests.md`.

## 2. Funding
- **Current text:** "Funding sources will be declared prior to submission."
- **Why it blocks:** placeholder wording; BMC requires a funding statement (including "no funding" if applicable).
- **Author confirmation needed:** list funders + grant numbers + funder role, or confirm no funding.
- **Proposed final text (template):** "This work received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors." *(Use only if true.)* — or a specific funding statement.
- **Who fills:** **Author input required.** Placeholder **removed from manuscript body** (v10.6).

## 3. CRediT author contributions
- **Current text:** "Author contributions will be finalized prior to submission."
- **Why it blocks:** placeholder; BMC requires per-author CRediT roles.
- **Author confirmation needed:** assign the 14 CRediT roles across the five authors; confirm "all authors read and approved the final manuscript."
- **Proposed final text:** complete the table in `author_contributions.md` (authors listed; roles to be filled).
- **Who fills:** **Author input required.** Placeholder **removed from manuscript body** (v10.6); template in `author_contributions.md`.

## 4. Acknowledgements
- **Current text:** "Acknowledgements will be added prior to submission."
- **Why it blocks:** placeholder; reads as unfinished.
- **Author confirmation needed:** list any non-author contributors / computing resources, or confirm "Not applicable."
- **Proposed final text:** "Not applicable." *(Use only if the authors confirm there is nothing to acknowledge.)*
- **Who fills:** **Author input required** (a one-line confirmation suffices). Placeholder **removed from manuscript body** (v10.6).

## 5. Paul Tan affiliation
- **Current text (byline, manuscript):** "Paul Tan³ … ³Harvard University Extension School."
- **Why it blocks:** earlier versions had "Affiliation to be finalized"; that is now resolved in text, but the affiliation has not been *confirmed* as accurate/final.
- **Author confirmation needed:** Paul Tan (and the corresponding authors) confirm that "Harvard University Extension School" is the correct, final affiliation and complete department/address as the journal requires.
- **Proposed final text:** keep as stated once confirmed; otherwise correct to the accurate affiliation.
- **Who fills:** **Author input required** (confirmation only). No body placeholder remains; nothing invented.

## 6. ORCIDs
- **Current state:** No ORCIDs in the manuscript byline; `CITATION.cff` carries commented "[TO BE SUPPLIED]" ORCID placeholders for all five authors.
- **Why it blocks:** BMC strongly encourages/links ORCIDs for all authors; missing ORCIDs are an open submission item.
- **Author confirmation needed:** each author supplies their ORCID iD (or confirms they will register one).
- **Proposed final text:** populate the `orcid:` fields in `CITATION.cff` and the author list at submission.
- **Who fills:** **Author input required.** No invented IDs.

## 7. Zenodo DOI / archived release
- **Current text (manuscript v10.5):** "…available in the project GitHub repository. The repository is structured for versioned archival, and a frozen DOI-citable release will accompany the final submission." (No DOI claimed — already cautious.)
- **Why it blocks (mildly):** the data-availability statement should give the actual repository URL and must not imply a DOI exists before one is minted.
- **Author/release confirmation needed:** deposit the v1.0.0 release to Zenodo and mint the DOI; then backfill it.
- **Proposed final text (applied):** "…available in the project GitHub repository, which is publicly available at https://github.com/iry496/leakage-aware-smo-svm-omics. A frozen DOI-citable archive will be deposited before final journal submission."
- **Who fills:** **Claude filled** the cautious wording + URL in the manuscript body (v10.6) and in `data_availability_statement.md`. **Author/release action required** to mint the DOI and insert it (in manuscript, `data_availability_statement.md`, `CITATION.cff`, release notes).

## 8. Cover letter declarations paragraph
- **Current text (v1):** "…All authors have approved the submission. Author contributions, competing-interest declarations, funding statements, and acknowledgements are being finalized and will be provided with the submission; ORCIDs will be supplied for all authors."
- **Why it blocks:** advertises that core declarations are unfinished — not a submission-ready letter.
- **Author confirmation needed:** none to clean the prose; the underlying declarations (items 1–4, 6, 7) remain unresolved, so the letter as a whole is not yet submission-ready.
- **Proposed final text (applied):** removed the "being finalized…" sentence; the letter now ends the originality paragraph at "All authors have approved the submission." A red **"DRAFT — NOT SUBMISSION-READY"** banner was added at the top listing the items to resolve before sending.
- **Who fills:** **Claude filled** the prose cleanup and the not-ready flag. The banner must be removed by the authors only after items 1–4, 6, 7 are resolved. **Cover letter remains NOT submission-ready.**

## 9. Data availability statement
- **Current text (v1, `data_availability_statement.md`):** Zenodo line read "DOI to be minted on the v1.0.0 release — [PLACEHOLDER: insert Zenodo DOI once published]."
- **Why it blocks:** placeholder DOI line; should use cautious, non-asserting wording until a DOI exists.
- **Author confirmation needed:** none for the wording; DOI minting is a release action (item 7).
- **Proposed final text (applied):** "The repository is publicly available at https://github.com/iry496/leakage-aware-smo-svm-omics. A frozen DOI-citable archive will be deposited before final journal submission." plus an explicit "[AUTHOR/RELEASE ACTION: mint Zenodo DOI and backfill]" note.
- **Who fills:** **Claude filled** the cautious wording. DOI insertion is an author/release action.

---

## Summary

**Claude completed safely (no invention):**
- Item 7 / 9: cautious data-availability wording + GitHub URL in the manuscript body (v10.6) and in `data_availability_statement.md`; no DOI asserted.
- Item 8: cleaned the cover letter prose and added a clear NOT-SUBMISSION-READY banner.
- Manuscript body: removed the four "will be declared/finalized/added prior to submission" placeholder subsections (Competing interests, Funding, Authors' contributions, Acknowledgements). Ethics, Consent, and Availability subsections remain and are complete.

**Unresolved — author input required (submission-blocking):**
1. Competing interests — declare or confirm none (all 5 authors).
2. Funding — declare funders/grants or confirm none.
3. CRediT author contributions — assign roles; confirm all-authors-approved.
4. Acknowledgements — supply or confirm "Not applicable."
5. Paul Tan affiliation — confirm "Harvard University Extension School" is final/complete.
6. ORCIDs — supply for all five authors.
7. Zenodo DOI — mint at release and backfill across files.
8. Cover letter — remove the DRAFT banner only after the above are resolved; confirm corresponding-author email and manuscript ID.
- Plus the two Step-4 citation fixes (GSE41998 → Horak 2013; OmicSelector preprint/software), pending approval.

**Confirmation:** No scientific result values were changed (numeric tokens identical between v10.5 and v10.6). No new models, feature selectors, elastic-net results, or BioTrust/founder/investor language were added.
