# GSE41998 External-Validation Design Memo

_Predeclared design memo. The protocol was subsequently implemented in
`scripts/14_external_validation_gse41998.R` and executed without changing its endpoint mapping,
feature mapping, scaling hierarchy, or go/no-go safeguards._

Discovery cohort: GSE25055 (Affymetrix HG-U133A, GPL96). Existing same-study-family validation:
GSE25065 (GPL96). Proposed new external cohort: **GSE41998** (different platform, GPL571 /
HG-U133A 2.0). This is the cross-platform test the manuscript currently lists as future work.

---

## 1. Dataset rationale
GSE25065 is a *same-platform, same-study-family* (Hatzis lineage) validation, so it tests
transportability only weakly. GSE41998 is an **independent, different-platform** neoadjuvant
breast-cancer cohort with pCR/RD-type response labels. Adding it would:
- provide the first genuinely **cross-platform** external test of the frozen guarded model;
- separate "transportability within a study family" from "transportability across platforms";
- strengthen the audit's external-validity claim without changing the within-cohort findings.
It remains an **audit/transportability** exercise, not biomarker discovery.

## 2. Label handling
- Include **only unambiguous pCR vs no-pCR (RD)** samples, coded to the same `pCR`/`RD` factor
  used for GSE25055/GSE25065 (positive class = `pCR`).
- **Exclude** samples coded `"0"`, blank, `NA`, `"NaN"`, `"N/A"`, or any unclear/ambiguous
  response value; exclude non-response arms or samples without an interpretable endpoint.
- Record, before any modeling, the counts kept/excluded and the resulting class balance; freeze
  the label mapping in code before predictions are generated.
- No relabeling or threshold-driven label redefinition after seeing performance.
- **Expected evaluable N (predeclared, from the dataset audit `tables/table1_dataset_audit_filled.csv`):**
  of 279 GSE41998 samples, the unambiguous endpoint is pCR = `Yes` (69) vs `No` (184); the **20 cases
  coded `"0"` and the 6 missing/unclear cases (26 total) are excluded as not-evaluable**, giving an
  expected evaluable **N = 253**. These exclusions stand **unless a documented GEO/publication data
  dictionary proves the `"0"` code has a definite pCR/RD meaning**. The script must report the realized
  kept/excluded counts and stop for review if they diverge materially from this predeclared N.

## 3. Platform mapping (GPL96 → GPL571)
- **Primary: exact probe-ID intersection.** Use only Affymetrix probe set IDs present in *both*
  platforms; project the frozen discovery model onto that shared probe space. Report how many of
  the model's features (the frozen discovery feature set) are recoverable in the intersection and
  how many are missing; if too few of the model's features survive (see go/no-go), do not force it.
- **Avoid gene-symbol collapse in the primary analysis** (collapsing many-to-one probe→gene
  introduces aggregation choices that can leak or distort). Use gene-symbol collapse **only as a
  pre-registered sensitivity analysis**, with a training-derived probe→gene map (max-variance or
  mean per gene fixed on discovery) and no use of validation data to choose the mapping.
- **Missing frozen probes:** the frozen feature set is the 100 discovery probes in
  `results/external_validation_gse25065/final_gse25055_selected_features.csv`. Report **by probe ID
  exactly which of these 100 are absent** from the GPL96∩GPL571 intersection, and **stop for review**
  rather than silently dropping them — **unless the predeclared fallback in §8 applies** (≥ ~70% of the
  100 frozen probes recovered ⇒ proceed on the recovered subset, documented; < ~70% ⇒ report-only, or
  the pre-registered gene-symbol sensitivity only).
- The mapping (intersection list or probe→gene rule) is **fixed on discovery before** GSE41998 is
  touched for prediction.

## 4. Scaling policy
- **Primary: discovery-derived scaling.** Apply the center/scale parameters estimated on the full
  GSE25055 discovery cohort (selected features only) unchanged to GSE41998 — identical to the
  GSE25065 protocol. No refitting on the external cohort.
- **Sensitivity: within-cohort, label-blind z-score.** Standardize GSE41998 features using that
  cohort's own mean/SD computed **without using the labels** (unsupervised), to gauge sensitivity
  to platform-scale shifts. This is explicitly an **unsupervised adaptation sensitivity analysis — NOT
  the frozen primary external validation** — and is reported and labeled separately so it is never
  mistaken for the strict-ML discovery-frozen result.
- Both options are specified in advance; neither is chosen post-hoc based on which gives a better
  number.

## 5. Leakage safeguards
- **No joint normalization** of discovery and validation (no pooled RMA/quantile normalization
  across cohorts); each cohort is processed independently and the model is frozen before contact.
- **No ComBat / batch correction across discovery and validation.** Cross-cohort batch correction
  that uses validation samples (especially with labels) is prohibited.
- **The full frozen model is fixed on discovery before any GSE41998 label is evaluated.** The frozen
  quantities are, with their committed discovery values: the **final feature set** (the 100 probes
  above; t-test top-K = 100 selection rule); the **SVM cost C = 0.25** (chosen by discovery-only
  guarded CV over {0.25, 1, 4}); the **class weights** (inverse class frequency, w_k = N/(K·n_k), from
  the discovery training labels); the **probability procedure** (e1071/libsvm built-in Platt-style
  estimates, fit on discovery); the **decision threshold** 0.5 on P(pCR); and the **center/scale
  parameters** (discovery-derived). **None of these — K, C, class weights, probability/scaling, mapping,
  threshold, or feature selection — is retuned or re-selected after seeing GSE41998 performance.**
- GSE41998 labels are used **only** for final evaluation, never for selection, scaling, mapping,
  or tuning.
- The mapping and scaling decisions are committed (script + memo) before predictions are produced;
  a single frozen prediction pass is run.

## 6. Expected outputs (if approved)
- `scripts/14_external_validation_gse41998.R` (canonical, frozen-model projection; mirrors
  `scripts/05_external_validation_gse25065.R`).
- `results/external_validation_gse41998/gse41998_external_metrics.csv` — AUROC, PR-AUC, balanced
  accuracy, MCC, sensitivity, specificity, confusion counts; primary (discovery-scaling) and the
  z-score sensitivity variant.
- `results/external_validation_gse41998/probe_mapping_report.csv` — intersection size, model
  features recovered/missing, mapping rule used.
- `results/external_validation_gse41998/external_validation_gse41998_notes.md` — provenance,
  label counts, safeguards checklist, cautious interpretation.
- `figures/external_validation_gse41998/` — discovery → GSE25065 → GSE41998 transportability
  comparison.
- An added row/panel to the evidence-audit artifact/dashboard (cross-platform transportability).

## 7. Risks
- **Probe-intersection shrinkage:** too few of the frozen model's features survive GPL96∩GPL571,
  making the projection underpowered or ill-defined.
- **Platform/scale shift:** different platform dynamic range can degrade the discovery-scaled
  projection (the z-score sensitivity variant is the hedge).
- **Label ambiguity:** GSE41998 response coding may differ; mis-mapping `"0"`/unclear values would
  corrupt the endpoint. Mitigated by strict include/exclude rules (§2).
- **Class imbalance / small usable n** after exclusions, widening uncertainty (report CIs).
- **Over-interpretation:** a low external number is a transportability limit, not evidence of
  leakage; framing must stay diagnostic/within-scope.

## 8. Go / no-go criteria
Proceed to run **only if all** hold (decided before generating predictions):
- **Labels:** ≥ ~80 usable samples with unambiguous pCR/RD after exclusions, and both classes
  present with pCR prevalence broadly comparable to discovery/GSE25065 (not, e.g., <5% or >60%).
- **Mapping:** the exact probe intersection recovers a clear majority (target ≥ ~70%) of the
  frozen model's features; if below that, do **not** force the primary analysis — either report it
  transparently as underpowered or fall back to the pre-registered gene-symbol sensitivity only.
- **Safeguards:** the frozen-model checklist (§5) is satisfied and committed before any GSE41998
  prediction.
- **No-go / report-only:** if labels are too few/ambiguous, or the intersection is too small, or
  the design cannot be executed without joint normalization/ComBat/retuning, do **not** run the
  primary projection; document the blocker and keep GSE41998 as stated future work.

## 9. Interpretation (predeclared)
- GSE41998 is a **cross-platform transportability sensitivity cohort**, not a confirmatory cohort.
- It is **not the main proof of leakage.** The leakage finding rests on the within-cohort
  naive-vs-guarded gap, the 1000-permutation negative control, and the 30-seed repeated CV;
  GSE41998 speaks only to cross-platform transportability.
- **Poor performance on GSE41998 does not invalidate the audit.** A low cross-platform number would
  demonstrate that **leakage control and transportability are distinct** questions, consistent with the
  paper's framing; it is **not evidence that leakage caused the external decline**.
- Required wording when reported: "cross-platform transportability sensitivity," "same-study-family
  validation" (for GSE25065), "supports the interpretation," "not biomarker discovery."

---

_This memo defines the protocol and its guardrails. The canonical script
`scripts/14_external_validation_gse41998.R` implements exactly this design (exact probe intersection,
discovery-frozen model, predeclared label rule, primary discovery-scaling + clearly separated
unsupervised z-score sensitivity, runtime go/no-go that stops report-only on failure) and writes no
raw expression. The executed run retained 253 explicit Yes/No cases and excluded 20 literal zeros
plus 6 missing values; the live GEO Series Matrix verification was repeated on 2026-08-07._
