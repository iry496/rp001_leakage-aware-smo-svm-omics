# Dataset Audit Report (v1)

Task 1 of the leakage-aware SMO/SVM omics study. This report records the metadata audit, endpoint label counts, sample-overlap assessment, and cohort-role recommendation for the five candidate GEO cohorts. No expression matrices or raw CEL files were downloaded; only GEO sample-characteristics metadata was parsed. No modeling was performed.

## Cohorts and verified label counts

Counts were tallied from each series' per-sample characteristics field (the exact field name used is recorded as the label source).

GSE25055 (discovery): platform GPL96 (HG-U133A), total N = 310. Label field `pathologic_response_pcr_rd`: pCR = 57, RD = 249, plus 4 samples coded NA (excluded). Status: verified from sample metadata. Source: Hatzis et al., JAMA 2011 (PMID 21558518).

GSE25065 (external validation 1): platform GPL96, total N = 198. Label field `pathologic_response_pcr_rd`: pCR = 42, RD = 140, plus 16 samples coded NA (excluded). Status: verified from sample metadata. Source: Hatzis et al., JAMA 2011.

GSE20194 (sensitivity; demoted from validation): platform GPL96, total N = 278. Label field `pcr_vs_rd`: pCR = 56, RD = 222, all samples labeled. Status: verified from sample metadata. Source: Shi et al., Nat Biotechnol 2010 / Popovici et al., Breast Cancer Res 2010.

GSE41998 (external validation 2): platform GPL571 (HG-U133A 2.0), total N = 279. Label field `pcr`: pCR (Yes) = 69, RD (No) = 184, used as the binary labels on 253 evaluable cases. A further 20 samples coded "0" and 6 with the field missing (26 total) are excluded from the binary endpoint. Status: verified exclusion rule from the live GEO Series Matrix on 2026-08-07; the exact clinical definition of the depositor's "0" code is not documented. Source: Horak et al., Clin Cancer Res 2013 (PMID 23340299).

GSE20271 (optional sensitivity): platform GPL96, total N = 178. Label field `pcr or rd`: pCR = 26, RD = 152, all samples labeled. Status: verified from sample metadata. Source: Tabchy et al., Clin Cancer Res 2010 (PMID 20829329).

## Endpoint harmonization

All five cohorts report a binary pathologic-complete-response versus residual-disease endpoint, so harmonization is feasible in principle. However the label is stored under a different field name in every cohort (`pathologic_response_pcr_rd`, `pcr_vs_rd`, `pcr or rd`, `pcr`), and value encodings differ (pCR/RD/NA versus Yes/No versus Yes/No/0). A single harmonization mapping must be defined and applied per cohort before any modeling.

## GSE41998 pCR coding verification

The live GSE41998 Series Matrix was re-read through GEOquery on 2026-08-07. The `pcr` field contains Yes (69), No (184), literal "0" (20), and missing (6). The explicit Yes and No cases are used as pCR and RD respectively. All 20 zeros remain zeros in the parallel `pcrrcb1` field, confirming that they are a distinct third category rather than a numeric recoding of No. The 20 zeros and 6 missing cases are therefore excluded from the binary endpoint. The exact clinical meaning of "0" remains undocumented in GEO and the accessible publication record.

## Open issue 2: MDACC FNA overlap risk

GSE20194 and GSE20271 are both MDACC fine-needle-aspiration cohorts drawn from overlapping neoadjuvant accrual programs, sharing investigators with the Hatzis discovery/validation cohorts (GSE25055 / GSE25065, SuperSeries GSE25066). Patient-level overlap among these series is plausible and well documented in the literature. Treating any two of them as independent would silently leak information into external validation. GSE25055 and GSE25065 are non-overlapping with each other by design (the original discovery/validation split), so that pairing is safe.

## Open issue 3: sample-level de-duplication required

Series-level metadata cannot confirm independence. Before GSE20194 or GSE20271 is used even as a sensitivity cohort, a sample-level de-duplication must be performed against the discovery and validation cohorts using clinical identifiers, hybridization/biopsy dates, and GSM-level matching. This step has not been done.

## GSE41998 label-resolution note

Goal: interpret the `pcr` sample-characteristics field (values Yes / No / 0 / missing) without downloading raw CEL files or modeling. Only GEO sample metadata and the publication abstract (Horak et al., Clin Cancer Res 2013, PMID 23340299) were consulted.

Distribution across the 279 samples: pcr = Yes 69, No 184, "0" 20, missing 6.

Cross-tabulation findings. The `pcr` field and the parallel `pcrrcb1` field agree exactly on Yes (69-69) and on "0" (20-20); they differ only within the No group, where pcrrcb1 reclassifies 17 of the 184 No cases as RCB-I responders (pcrrcb1 is the broader pCR-or-RCB-I endpoint). The "0" value is therefore a distinct third category, not a numeric 0/1 recoding of No. By treatment arm, the "0" cases span all arms (6 ixabepilone, 6 paclitaxel, 8 non-randomized "none"), and the 6 missing-pcr cases are all in the non-randomized "none" arm. The separate `ac response` field (clinical response to AC induction: complete/partial/stable/progressive/unable-to-determine) is an independent axis: "0" pcr cases include 4 complete, 12 partial, 2 stable, 1 progressive, and 1 unable-to-determine AC responses, so "0" does not map onto any single clinical-response category. (This corrects the earlier v1 note, which conflated AC-phase clinical response with pathologic pCR.)

Interpretation. The data establish that "0" is not equivalent to the explicit No category: it is preserved separately in both `pcr` and `pcrrcb1` and occurs across treatment arms and AC-response categories. The conservative and reproducible analysis rule is therefore to exclude "0" from the pCR-versus-RD endpoint rather than infer a clinical label. Neither the GEO series record nor the accessible publication record provides a coding key that defines the clinical meaning of "0".

Decision. pCR_N = 69 and RD_N = 184 are recorded from the explicit Yes/No labels (253 evaluable cases). The 20 "0" and 6 missing cases (26 total) are excluded from the binary counts. Count_status for GSE41998 is "verified exclusion rule; exact zero-code meaning undocumented." This uncertainty does not affect the analyzed 253 explicit Yes/No cases.

## Cohort-role recommendation

Discovery: GSE25055. External validation 1: GSE25065 (same platform; minimal harmonization). External validation 2: GSE41998 (independent trial, on platform GPL571 requiring cross-platform harmonization; analyze on the 253 evaluable Yes/No cases). Optional sensitivity: GSE20271 (only after de-duplication). GSE20194 is demoted from external validation to sensitivity-only because of its overlap risk with the Hatzis cohorts.

## Constraints honored

No raw CEL files or expression matrices downloaded; only GEO sample-characteristics metadata and the publication abstract parsed. No modeling performed. protocol_v0_1.md not edited.
