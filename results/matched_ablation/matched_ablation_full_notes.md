# Matched one-factor ablation - full

- GSE25055: n=306; repeated nested CV: 30 x 5 outer folds, 5 inner folds.
- The only unmatched component is supervised feature-selection placement.
- All outer test-fold keys passed exact arm-identity checks.
- Final-model RNG seeds were identical across arms within every repeat and outer fold.
- Fast Welch selector passed exact-order equivalence against the canonical selector.

## Paired split-level results
- AUROC: leaky mean 0.7870; guarded mean 0.7256; mean delta 0.0615 (empirical split interval 0.0084 to 0.1201; positive in 30/30 repeats).
- PR-AUC: leaky mean 0.4237; guarded mean 0.3733; mean delta 0.0505 (empirical split interval -0.0227 to 0.1310; positive in 24/30 repeats).
- Split repeats reuse the same patients and are correlated; the empirical interval is a robustness distribution, not an independence-based confidence interval.
- Patient-bootstrap intervals use repeat-averaged cross-fitted scores and are conditional on the realized CV fits; the workflow is not refit within bootstrap samples.

## Interpretation
This matched contrast isolates supervised feature-selection placement under the stated pipeline contract. It replaces the older unmatched complete-workflow contrast as the leakage-effect estimate; the older comparison remains useful only as a complete-workflow audit.

_Runtime 1.7 minutes._
