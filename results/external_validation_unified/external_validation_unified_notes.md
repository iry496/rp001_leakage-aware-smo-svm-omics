# Unified frozen-model external validation

- One GSE25055 model was fit once and reused for all projections (top-100; cost=4).
- Both primary external projections use the identical discovery-derived scaler and explicit P(pCR) >= 0.50 rule.
- Package-class versus explicit-threshold discordance: 0 predictions.
- Frozen feature recovery: GSE25065 100/100; GSE41998 100/100.

## Primary external results
- GSE25065 (n=182; pCR prevalence 0.231): AUROC 0.5633, PR-AUC 0.2986, Brier 0.1970, calibration intercept -0.927, slope 0.148; sensitivity 0.119, specificity 0.943.
- GSE41998 primary (n=253; pCR prevalence 0.273): AUROC 0.6922, PR-AUC 0.4634, Brier 0.2668, calibration intercept -1.229, slope 0.617; sensitivity 0.797, specificity 0.500.
- GSE41998 label-blind z-score sensitivity: AUROC 0.7067, PR-AUC 0.4674, Brier 0.1873; sensitivity 0.101, specificity 0.962.

## Interpretation and uncertainty
Threshold-independent discrimination is primary. PR-AUC must be interpreted against each cohort's pCR prevalence.
The same 0.5 operating point behaves very differently across cohorts and scaling variants; this is evidence of transport/calibration instability, not a reason to retune on external outcomes.
Bootstrap intervals are conditional on the single frozen model and stored external predictions. They quantify external sample variation, not discovery-model refitting uncertainty.
GSE41998 combines platform, treatment, population, and endpoint-context transport; its result must not be attributed to platform alone.
Legacy rounded-result regression checks: 5/5 passed within tolerance.

_Runtime 2.1 minutes._
