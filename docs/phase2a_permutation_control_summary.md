# Phase 2A Completion Summary — Permutation Control

_Status snapshot for the BioData Mining methodology manuscript_
**"An Audit of Data Leakage and Feature Stability in Transcriptomic Biomarker Classification Using Nested SMO/SVM and External Validation."**

This summary reflects only material already merged into `main` (PR #24). It introduces no new analyses, manuscript edits, or code beyond what is committed, and changes no existing result values. All numbers below are read directly from the committed B=1000 outputs.

---

## 1. What Phase 2A completed

Phase 2A added a **label-shuffle negative control** (permutation test) for the leakage audit on the discovery cohort GSE25055. The control re-runs **both** pipelines — the leaky baseline (global feature selection before cross-validation) and the guarded nested SMO/SVM (feature selection inside training folds only) — on randomly permuted pCR/RD labels, building a null distribution for each performance metric and an empirical p-value for the observed (real-label) result.

The purpose is diagnostic: to show how each pipeline behaves when there is, by construction, no real signal. A guarded pipeline should collapse toward chance under shuffled labels; a leaky pipeline should not. This directly addresses Phase 1 limitation #1 (the leakage gap was positive but not statistically significant in a single cohort), by characterising the leakage mechanism rather than relying on the single-cohort gap.

## 2. PR merged (#24)

| PR | Branch | Merge commit | Summary |
|----|--------|--------------|---------|
| #24 | `analysis/permutation-control` | `39c159b` | B=1000 label-shuffle permutation control (GSE25055 only): canonical R script, optional Python cross-check, null-distribution table, per-statistic p-values, notes, and figure. |

The branch remains on the remote (not deleted). The PR contained exactly eight files and no manuscript, external-cohort, raw, processed, `.rds`/`.RData`, GEO, rendered-HTML, or cache content.

## 3. Analysis design

- **Cohort:** GSE25055 only (306 usable samples; 57 pCR, 249 RD; pCR prevalence 0.1863). External cohorts (GSE25065, GSE41998, GSE20194, GSE20271) are explicitly **not** loaded — enforced by a scope-guard in the script.
- **Permutations:** identity (permutation 0) + **1000** shuffled-label replicates. Only the sample labels are permuted; the expression matrix is never altered.
- **Leaky arm:** univariate Welch t-test, top-K = 100 probes selected on the **full** data **before** cross-validation, then 5×5 repeated stratified CV.
- **Guarded arm:** nested CV (5 outer × 5 inner folds); feature selection (same Welch t-test, top-100) and cost tuning (grid 0.25/1/4) performed **inside training folds only**.
- **Model / metrics:** linear SMO/SVM (libsvm via e1071, class-weighted, probability outputs); AUROC (pROC), PR-AUC (PRROC `auc.integral`), plus balanced accuracy, MCC, sensitivity, specificity.
- **Folds** are regenerated (stratified on the shuffled labels) and feature selection re-run for every permutation, so no information from the real labels leaks into the null.
- **Seeds:** master seed 20260620; per-permutation label shuffles use offset +100000. Identity must reproduce the committed pilot values.

## 4. Integrity checks (identity + selector)

Two guardrails were verified on every run and are recorded in the committed notes:

- **Identity reproduction (PASS):** with unpermuted labels the script reproduces the pilot exactly — leaky AUROC **0.7705** (ref 0.7705) and guarded AUROC **0.7265** (ref 0.7265), both within the 0.01 tolerance. This confirms the permutation harness is a faithful wrapper around the original pipelines.
- **Selector validation (PASS):** the optimized vectorized Welch selector was compared against the original `t.test`-based selector on the identity case plus five shuffles. Result: **100/100 probe overlap and identical ranking (exact_order = TRUE)** in all six cases. The speed optimization is therefore scientifically inert.

## 5. Key results — AUROC null distributions (B = 1000)

| Statistic | Observed (real labels) | Null mean | Null 2.5% | Null 97.5% | Fraction of null > chance |
|-----------|-----------------------:|----------:|----------:|-----------:|--------------------------:|
| Leaky AUROC | 0.7705 | 0.8783 | 0.7613 | 0.9564 | 1.000 |
| Guarded AUROC | 0.7265 | 0.5401 | 0.4783 | 0.6185 | 0.871 |
| Leakage gap (leaky − guarded) | 0.0440 | 0.3382 | 0.2090 | 0.4454 | 1.000 (gap > 0) |

Under shuffled labels the **guarded** pipeline sits close to chance (null mean 0.54), while the **leaky** pipeline remains far above chance (null mean 0.88, with every one of the 1000 shuffles exceeding 0.5). The null leakage gap is large and entirely positive.

## 6. PR-AUC null distributions (reference = prevalence 0.1863)

| Statistic | Observed | Null mean | Null 2.5% | Null 97.5% | Fraction above prevalence |
|-----------|---------:|----------:|----------:|-----------:|--------------------------:|
| Leaky PR-AUC | 0.4020 | 0.6637 | 0.4131 | 0.8606 | 1.000 |
| Guarded PR-AUC | 0.3656 | 0.1917 | 0.1465 | 0.2605 | 0.489 |

The guarded PR-AUC null centres on the cohort prevalence (0.192 vs 0.1863), with roughly half the replicates above and half below — the expected, well-calibrated behaviour for a no-signal model. The leaky PR-AUC null inflates far above prevalence.

## 7. Empirical p-values and interpretation

`p_obs_ge_null` is the one-sided permutation p-value (fraction of null replicates at or above the observed value):

| Statistic | Observed | p(obs ≥ null) |
|-----------|---------:|--------------:|
| Leaky AUROC | 0.7705 | 0.967 |
| Guarded AUROC | 0.7265 | **0.001** |
| Leakage gap (AUROC) | 0.0440 | 1.000 |
| Leaky PR-AUC | 0.4020 | 0.979 |
| Guarded PR-AUC | 0.3656 | **0.001** |
| Leakage gap (PR-AUC) | 0.0363 | 1.000 |

Reading these together:

- The **guarded** pipeline's real-label performance significantly exceeds its own shuffled-label null (p = 0.001 for both AUROC and PR-AUC). This supports the interpretation that the guarded model captures genuine label-associated signal, and that the guarded null is well-calibrated near chance.
- The **leaky** pipeline reports high performance even from random labels (AUROC ≈ 0.88), and its real-label score falls **inside / below** its own noise distribution (p ≈ 0.97). This supports the interpretation that feature-selection leakage can inflate apparent performance independent of any real signal. It is a diagnostic of the leakage mechanism, **not** evidence of biological signal.
- The observed real-data leakage gap (0.044) sits well below the null leakage gap (mean 0.338). The permutation control therefore demonstrates the leakage **mechanism** clearly (the two nulls diverge sharply) even though the two arms land close together on this particular real dataset — consistent with the Phase 1 finding that the single-cohort gap was not statistically significant.

## 8. Optimization and runtime

To make B = 1000 feasible, the permutation driver was optimized without changing the science:

- **Vectorized Welch selector** (matrix-algebra t-statistic) replacing per-probe `t.test`, validated as bit-identical in selection (Section 4).
- **Parallel execution** across permutations (`parallel::mclapply`, 4 workers), with serial and original-selector fallbacks retained via toggles.
- **RDS caching** of the processed matrix/labels (git-ignored), so the matrix is loaded once.

Runtime improved from **361.8 s per permutation** (serial, original selector) to **6.2 s compute / 1.6 s wall per permutation** (4 workers). The full B = 1000 run completed in **≈ 26.6 minutes**. The optimized B=200 checkpoint produced null distributions bit-identical to the original slow smoke run (max absolute difference 0.0 across all metric cells).

## 9. Reproducibility and scope guards

- **Seeds and fold assignments** are fixed (master 20260620; permutation offset +100000); RNG is reset inside the run functions for determinism.
- **Scope guard** in the script prevents loading any external cohort.
- **Data hygiene:** `.gitignore` excludes `processed_data/`, `raw_data/`, and R binary caches (`*.rds`, `*.RData`, `*.CEL`, `*.gz`); none were committed.
- **Canonical results are from R.** The Python script (`scripts/checks/08_permutation_control_smoke_check.py`) is a gated, directional cross-check only (different SVM solver and fold RNG); it requires a local matrix export and never fetches data or fabricates results.

## 10. Files added in Phase 2A (on `main`)

| File | Role |
|------|------|
| `scripts/08_permutation_control.R` | Canonical permutation-control script (both arms, identity, selector validation, parallel, CLI `RUN_TAG B_PERM`). |
| `scripts/checks/08_permutation_control_smoke_check.py` | Optional directional Python cross-check (gated on a local matrix export). |
| `results/permutation/permutation_b1000_null_distributions.csv` | Per-permutation metrics (identity + 1000 shuffles). |
| `results/permutation/permutation_b1000_pvalues.csv` | Per-statistic observed value, null summary, and empirical p-value. |
| `results/permutation/permutation_b1000_notes.md` | Human-readable run notes (design, integrity checks, null summary, runtime). |
| `figures/permutation_b1000_null.png` / `.pdf` | Null AUROC histograms (leaky vs guarded, with observed lines) and null leakage-gap histogram. |
| `.gitignore` | Updated to exclude raw/processed data and R binary caches. |

## 11. Remaining limitations and what is still needed

1. **Single cohort.** The permutation control is on GSE25055 only; it characterises the leakage mechanism within one dataset, not across cohorts.
2. **Single selector / single classifier.** Univariate Welch t-test top-K and linear SMO/SVM only; the null behaviour of other selectors/classifiers is untested here.
3. **Real-data leakage gap remains modest.** The control shows the mechanism, but the observed gap on real labels is small and (from Phase 1) not statistically significant in this cohort. Repeated nested CV (#12) is still needed to characterise the gap and feature-stability distributions across seeds.
4. **R-execution caveat.** R is not runnable in the working/automation environment; the committed B=1000 outputs were generated by running the canonical R script in the user's RStudio/Terminal. Point estimates and identity reproduction are deterministic; the null is reproducible under the fixed seeds.
5. **Claim calibration not yet finalized.** The manuscript (#9 `manuscript/claim-calibration`) should be updated to incorporate this permutation evidence using cautious, diagnostic language — once repeated nested CV is also available.

## 12. Recommended next Phase 2 tasks

Per the revision roadmap (`docs/revision_roadmap_after_ai_review.md`):

- **#12 `analysis/repeated-nested-cv`** (P1, new modeling) — repeat nested CV across 20–30 seeds; report leakage-gap and feature-stability distributions. This pairs with the permutation control to address Phase 1 limitation #1.
- **#9 `manuscript/claim-calibration`** — fold the bootstrap CIs, permutation control, and repeated-CV results into calibrated claims.
- Later phases (unchanged): #13 `analysis/gse41998-external`, #14 `analysis/comparator-elasticnet`; then #15 `figure/evidence-audit-dashboard` and #16 `repo/reproducibility-release`; optional #17 `analysis/selector-and-k-sweep`, #18 `analysis/calibration-threshold`.

## 13. One-line takeaway

The B=1000 label-shuffle control shows the guarded nested pipeline behaving correctly under the null (near chance, real-label performance significantly above its null at p = 0.001) while the leaky pipeline produces strongly above-chance results from random labels — supporting the interpretation that feature-selection leakage can inflate apparent performance, as a diagnostic of the audited mechanism rather than evidence of biological signal.
