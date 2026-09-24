# Citation Verification — v1

**Manuscript checked against:** `manuscript/Leakage_Aware_SMO_SVM_Manuscript_v10_5_table_figure_readability.docx` (reference list paragraphs and in-text mentions).
**Date:** 2026-06-24
**Status:** Verification report only — **no manuscript edits made.** All wording recommendations below are proposals pending author approval. No unverified citations were added.

## Summary table

| # | Item | In manuscript? | Verified status | Action needed |
| --- | --- | --- | --- | --- |
| 1 | nestedcv / Lewis et al. 2023 | Cited (ref list + in-text) | Peer-reviewed + software | None — correct |
| 2 | OmicSelector | **Mentioned in text, NOT cited** | Preprint (bioRxiv) + software | **Add citation or adjust wording** |
| 3 | Horak et al. / GSE41998 | **Not cited at all** | Peer-reviewed | **Add as GSE41998 data-source citation** |
| 4 | Hatzis et al. / GSE25055/65 | Cited | Peer-reviewed | None — correct |
| 5 | TRIPOD+AI / Collins et al. 2024 | Cited | Peer-reviewed | None — correct |
| 6 | BioDiscML / BioDiscViz | **Not mentioned** | Peer-reviewed + software | Optional only (not currently referenced) |
| 7 | Nogueira et al. (feature stability) | Cited | Peer-reviewed | None (optional: add URL) |
| 8 | Kapoor & Narayanan (leakage) | Cited | Peer-reviewed | None — correct |
| 9 | Ambroise & McLachlan (FS leakage) | Cited | Peer-reviewed | None — correct |

---

## 1. nestedcv — Lewis et al. (2023)

- **Exact citation:** Lewis, M. J., Spiliopoulou, A., Goldmann, K., Pitzalis, C., McKeigue, P., & Barnes, M. R. (2023). nestedcv: an R package for fast implementation of nested cross-validation with embedded feature selection designed for transcriptomics and high-dimensional data. *Bioinformatics Advances*, 3(1), vbad048.
- **DOI / PMID:** DOI 10.1093/bioadv/vbad048 · PMID 37113250 · PMC10125905
- **Publication status:** Peer-reviewed journal article (Oxford University Press) describing a software tool (R package `nestedcv`, on CRAN).
- **Cited in manuscript?** Yes — reference list entry, plus in-text at the "Remaining methodological gap" and contribution paragraphs.
- **Wording change?** None. The manuscript's reference entry matches the verified record exactly (authors, year, journal, volume/issue, article id, DOI). The in-text framing ("nestedcv … provides nested cross-validation with embedded feature selection") is accurate.

## 2. OmicSelector — Stawiski et al. (2022)

- **Exact citation:** Stawiski, K., et al. (2022). OmicSelector: automatic feature selection and deep learning modeling for omic experiments. *bioRxiv* 2022.06.01.494299.
- **DOI:** 10.1101/2022.06.01.494299 (bioRxiv preprint). Software: R package / Docker application, GitHub `kstawiski/OmicSelector`.
- **Publication status:** **Preprint (bioRxiv, 2022) + software.** As of this check (June 2026) no separate peer-reviewed journal article describing the OmicSelector *method* was found; the tool has been *used* in later peer-reviewed studies, but the primary methodology reference remains the bioRxiv preprint plus the software repository.
- **Cited in manuscript?** **Mentioned but NOT cited.** OmicSelector appears in two places with a substantive capability claim — "automated omics environments such as OmicSelector provide leakage-resistant feature selection and modeling, including nested cross-validation and stability metrics" (Related work / methodological-gap paragraph) and "This work does not replace nestedcv, stabm, OmicSelector, or other existing tools …" (contribution paragraph). There is **no reference-list entry** for it.
- **Wording change? — ACTION REQUIRED.** Making a specific capability claim about a named tool without a citation is a desk-review risk. Two acceptable fixes (author to choose):
  - (a) Add a citation to the bioRxiv preprint + GitHub repository, and clearly mark it as a preprint/software (not a peer-reviewed methods paper); or
  - (b) Soften the claim to avoid asserting specific verified capabilities (e.g., "automated omics environments (e.g., OmicSelector)") if a citation is not desired.
  Recommended: option (a), with the preprint/software status stated, so the claim is attributable. Per instruction, no citation has been added yet.

## 3. Horak et al. (2013) — original publication for GSE41998

- **Exact citation:** Horak, C. E., Pusztai, L., Xing, G., Trifan, O. C., Saura, C., Tseng, L.-M., Chan, S., Welcher, R., & Liu, D. (2013). Biomarker analysis of neoadjuvant doxorubicin/cyclophosphamide followed by ixabepilone or paclitaxel in early-stage breast cancer. *Clinical Cancer Research*, 19(6), 1587–1595.
- **DOI / PMID:** DOI 10.1158/1078-0432.CCR-12-1359 · associated trial NCT00455533.
- **Publication status:** Peer-reviewed journal article. This is the study whose expression data are deposited as **GSE41998** (Affymetrix HG-U133A 2.0 / GPL571; AC followed by randomization to ixabepilone vs paclitaxel) — the cross-platform transportability sensitivity cohort.
- **Cited in manuscript?** **No.** GSE41998 is used and described 20+ times (Methods, Results, Table 11/12, Figure 7) but its **source publication is not cited** — there is no "Horak" entry anywhere in the manuscript.
- **Wording change? — ACTION REQUIRED.** For data provenance and desk-review compliance, the GSE41998 cohort should cite its original publication (Horak et al., 2013). Recommended: add the reference and cite it at first mention of GSE41998 in Methods (analogous to how GSE25055/65 cite Hatzis et al., 2011). Per instruction, the reference has not been added yet. (Note: the manuscript's treatment-regimen description — "AC-then-ixabepilone-or-paclitaxel" — already matches this source, so only the citation is missing, not the description.)

## 4. Hatzis et al. (2011) — GSE25055 / GSE25065

- **Exact citation:** Hatzis, C., Pusztai, L., Valero, V., Booser, D. J., Esserman, L., Lluch, A., Vidaurre, T., Holmes, F., Souchon, E., Wang, H., Martin, M., Cotrina, J., Gomez, H., Hubbard, R., Chacón, J. I., Ferrer-Lozano, J., Dyer, R., Buxton, M., Gong, Y., & Symmans, W. F. (2011). A genomic predictor of response and survival following taxane-anthracycline chemotherapy for invasive breast cancer. *JAMA*, 305(18), 1873–1881.
- **DOI / PMID:** DOI 10.1001/jama.2011.593 · PMID 21558518
- **Publication status:** Peer-reviewed journal article. Source of the GSE25066 SuperSeries and its non-overlapping subseries GSE25055 (discovery) and GSE25065 (validation).
- **Cited in manuscript?** Yes — reference list entry; cohorts described as discovery / same-study-family validation.
- **Wording change?** None. Reference matches the verified record exactly.

## 5. TRIPOD+AI — Collins et al. (2024)

- **Exact citation:** Collins, G. S., Moons, K. G. M., Dhiman, P., Riley, R. D., Beam, A. L., Van Calster, B., Ghassemi, M., Liu, X., Reitsma, J. B., van Smeden, M., … & members of the TRIPOD+AI initiative (2024). TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ*, 385, e078378.
- **DOI / PMID:** DOI 10.1136/bmj-2023-078378 · PMC11025451 (published 16 April 2024)
- **Publication status:** Peer-reviewed reporting-guideline statement (BMJ).
- **Cited in manuscript?** Yes — reference list entry.
- **Wording change?** None. Reference matches the verified record (authors, year, journal, article id, DOI).

## 6. BioDiscML / BioDiscViz

- **Exact citation (BioDiscML):** Leclercq, M., Vittrant, B., Martin-Magniette, M. L., Scott Boyer, M. P., Perin, O., Bergeron, A., Fradet, Y., & Droit, A. (2019). Large-scale automatic feature selection for biomarker discovery in high-dimensional OMICs data. *Frontiers in Genetics*, 10, 452.
- **DOI:** 10.3389/fgene.2019.00452. Software: GitHub `mickaelleclercq/BioDiscML`. **BioDiscViz** is a companion visualization/consensus-signature tool reported later (2023).
- **Publication status:** Peer-reviewed journal article (Frontiers in Genetics) + software.
- **Cited in manuscript?** **No — not mentioned at all** (zero occurrences of "BioDiscML" or "BioDiscViz" in the text).
- **Wording change?** None required, because the tool is not referenced. Optional only: if the author wants to broaden the "existing automated omics tools" comparison (alongside nestedcv and OmicSelector), BioDiscML would be a legitimate, verifiable addition — but this is discretionary and not a desk-review blocker. No action taken.

## 7. Nogueira et al. (2018) — feature-selection stability

- **Exact citation:** Nogueira, S., Sechidis, K., & Brown, G. (2018). On the stability of feature selection algorithms. *Journal of Machine Learning Research*, 18(174), 1–54.
- **DOI / locator:** JMLR has no DOIs; canonical URL https://jmlr.org/papers/v18/17-514.html
- **Publication status:** Peer-reviewed journal article (JMLR). Source of the Nogueira stability index used throughout the manuscript.
- **Cited in manuscript?** Yes — reference list entry; the "Nogueira index" is referenced ~19 times.
- **Wording change?** None substantive. The reference is correct (authors, year, volume/issue 18(174), pages 1–54). Optional nicety: add the JMLR article URL, since the entry currently has no DOI/URL (JMLR articles have none, so this is acceptable as-is).

## 8. Kapoor & Narayanan (2023) — leakage / reproducibility

- **Exact citation:** Kapoor, S., & Narayanan, A. (2023). Leakage and the reproducibility crisis in machine-learning-based science. *Patterns*, 4(9), 100804.
- **DOI:** 10.1016/j.patter.2023.100804
- **Publication status:** Peer-reviewed journal article (Patterns, Cell Press).
- **Cited in manuscript?** Yes — reference list entry; cited in the leakage-framing paragraphs.
- **Wording change?** None. Reference matches the verified record exactly.

## 9. Ambroise & McLachlan (2002) — feature-selection leakage

- **Exact citation:** Ambroise, C., & McLachlan, G. J. (2002). Selection bias in gene extraction on the basis of microarray gene-expression data. *Proceedings of the National Academy of Sciences*, 99(10), 6562–6566.
- **DOI / PMID:** DOI 10.1073/pnas.102102699 · PMID 11983868 · PMC124442
- **Publication status:** Peer-reviewed journal article (PNAS).
- **Cited in manuscript?** Yes — reference list entry; cited in the leakage-framing paragraphs.
- **Wording change?** None. Reference matches the verified record exactly.

---

## Bottom line (for the next manuscript-edit step, pending approval)

Two items are genuine, verifiable gaps that a desk reviewer could flag:

1. **GSE41998 has no source citation.** Add Horak et al. (2013), *Clin Cancer Res*, 19(6):1587–1595, DOI 10.1158/1078-0432.CCR-12-1359, and cite it at the first GSE41998 mention (parallel to Hatzis et al. for GSE25055/65).
2. **OmicSelector is claimed but uncited.** Either add the bioRxiv preprint + software citation (marked as preprint/software, not peer-reviewed), or soften the capability claim.

Everything else (nestedcv, Hatzis, TRIPOD+AI, Nogueira, Kapoor & Narayanan, Ambroise & McLachlan) is correctly cited and verified. BioDiscML/BioDiscViz is not mentioned and needs no action unless the author chooses to add it. No citations have been added or changed in the manuscript; awaiting approval before any edit.
