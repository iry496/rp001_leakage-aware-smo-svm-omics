# GSE41998 external validation (GO) — frozen GSE25055 model

- Status: GO. Usable n=253 (pCR=69, RD=184; prevalence 0.273). Label column: 'characteristics_ch1.11'.
- Frozen feature set: top-100 (Welch t-test on full GSE25055); recovered in GSE41998: 100/100 (coverage 1.000).
- Frozen cost=4 (guarded discovery CV); threshold=0.50 on P(pCR); no GSE41998 tuning/selection.
- Exact probe-ID intersection; no gene-symbol collapse; no joint normalization; no ComBat.

## Results (frozen-model projection)
- PRIMARY (discovery-derived scaling): AUROC 0.6922, PR-AUC 0.4634, balanced acc 0.6486, MCC 0.2682, sens 0.7971, spec 0.5000.
- SENSITIVITY (within-cohort z-score): AUROC 0.7067, PR-AUC 0.4674.

## Interpretation
Cross-platform transportability of a frozen model; a generalization limit, not a leakage effect.
Within-cohort/diagnostic; not biomarker discovery. GSE41998 labels used only for final evaluation.

_Runtime 1.6 min. No raw expression written._
