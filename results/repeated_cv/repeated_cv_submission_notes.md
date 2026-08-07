# Repeated Nested CV - submission run (GSE25055 only)

- Seeds (30): 20260620, 20260621, 20260622, 20260623, 20260624, 20260625, 20260626, 20260627, 20260628, 20260629, 20260630, 20260631, 20260632, 20260633, 20260634, 20260635, 20260636, 20260637, 20260638, 20260639, 20260640, 20260641, 20260642, 20260643, 20260644, 20260645, 20260646, 20260647, 20260648, 20260649
- Matched seed-level comparison (NOT identical fold structure): leaky = 5x5 repeated CV,
  guarded = nested 5-outer x 5-inner; within each seed both arms draw stratified folds
  from the same master seed.
- Fast Welch selector: TRUE ; parallel across seeds: TRUE (4 workers); CV loops serial.
- Leaky arm: global t-test top-100 before CV. Guarded arm: t-test inside training folds only.
- Cost grid: 0.25, 1, 4. Top-K = 100.

## Anchor reproduction check (seed 20260620)
- leaky AUROC = 0.7830 (ref 0.7830)
- guarded AUROC = 0.7032 (ref 0.7032)
- Nogueira = 0.5128 (ref 0.5128)
- mean Jaccard = 0.3487 (ref 0.3487)
- stable-core count = 26 (ref 26)
- unstable-tail count = 105 (ref 105)
- Overall: PASS

## Selector validation (original vs fast Welch)
- identity: overlap 100/100, exact_order=TRUE.
- perm1: overlap 100/100, exact_order=TRUE.
- perm2: overlap 100/100, exact_order=TRUE.
- perm3: overlap 100/100, exact_order=TRUE.
- perm4: overlap 100/100, exact_order=TRUE.
- perm5: overlap 100/100, exact_order=TRUE.

## Per-seed leakage gap (AUROC)
- leaky AUROC: 0.7830, 0.7767, 0.7726, 0.7740, 0.7716, 0.7666, 0.7630, 0.7600, 0.7620, 0.7686, 0.7724, 0.7816, 0.7838, 0.7781, 0.7781, 0.7789, 0.7691, 0.7646, 0.7661, 0.7659, 0.7580, 0.7528, 0.7607, 0.7608, 0.7547, 0.7583, 0.7688, 0.7628, 0.7551, 0.7582
- guarded AUROC: 0.7032, 0.7532, 0.7452, 0.7356, 0.7824, 0.7406, 0.7161, 0.6870, 0.7123, 0.6895, 0.7239, 0.6939, 0.6915, 0.6998, 0.7391, 0.7567, 0.7000, 0.7406, 0.7403, 0.7118, 0.6984, 0.6612, 0.7289, 0.7220, 0.7663, 0.7189, 0.7289, 0.7623, 0.7502, 0.7060
- gap AUROC (leaky-guarded): 0.0799, 0.0235, 0.0274, 0.0384, -0.0108, 0.0260, 0.0469, 0.0731, 0.0497, 0.0791, 0.0485, 0.0877, 0.0923, 0.0783, 0.0390, 0.0222, 0.0691, 0.0239, 0.0258, 0.0541, 0.0596, 0.0916, 0.0319, 0.0388, -0.0116, 0.0394, 0.0399, 0.0006, 0.0049, 0.0522
- gap PR-AUC: 0.0781, -0.0167, 0.0559, -0.0162, -0.0464, 0.0179, 0.0486, 0.0598, 0.0786, 0.0578, 0.0646, 0.0950, 0.0975, 0.0490, 0.0370, 0.0215, 0.0460, 0.0034, -0.0144, 0.0697, 0.0754, 0.0680, 0.0599, -0.0519, -0.0665, 0.0631, 0.0356, -0.0341, -0.0142, 0.0284

## Per-seed feature stability
- Nogueira: 0.5128, 0.5550, 0.5560, 0.5078, 0.5751, 0.5610, 0.5671, 0.5429, 0.5480, 0.5440, 0.5590, 0.5078, 0.5118, 0.5349, 0.5480, 0.5510, 0.5219, 0.5349, 0.5399, 0.5339, 0.5540, 0.4847, 0.5208, 0.5590, 0.5429, 0.5500, 0.6002, 0.5560, 0.4746, 0.5269
- mean Jaccard: 0.3487, 0.3881, 0.3892, 0.3440, 0.4065, 0.3921, 0.3989, 0.3764, 0.3796, 0.3761, 0.3908, 0.3445, 0.3481, 0.3697, 0.3798, 0.3840, 0.3572, 0.3680, 0.3727, 0.3674, 0.3863, 0.3229, 0.3548, 0.3906, 0.3766, 0.3824, 0.4310, 0.3891, 0.3185, 0.3604
- stable-core count: 26, 25, 28, 24, 33, 30, 31, 30, 30, 30, 31, 24, 23, 29, 32, 29, 24, 26, 26, 25, 32, 24, 30, 29, 28, 32, 35, 32, 17, 27
- unstable-tail count: 105, 109, 99, 133, 91, 98, 91, 113, 104, 111, 94, 105, 112, 116, 101, 101, 105, 105, 105, 102, 98, 131, 114, 96, 101, 122, 87, 109, 125, 116

## Runtime
- mean 9.0 s/seed (compute, within worker); wall 3.5 s/seed; total 1.8 min.
- estimated full 20 seeds: 1.2 min (0.02 h) ; 30 seeds: 1.8 min (0.03 h) (at current wall throughput).

## Gap significance (Wilcoxon signed-rank, paired delta vs 0)
- delta AUROC: V = 458, p = 3.69e-06; positive in 28/30 seeds.
- delta PR-AUC: V = 390, p = 0.00124; positive in 22/30 seeds.
- The workflow contrasts are positive in 28/30 seeds for AUROC and 22/30 seeds for PR-AUC (Wilcoxon p = 3.69e-06 and 0.00124, respectively), but the 2.5-97.5% intervals include small negative values. These are seed/fold robustness summaries, not patient-level inference or estimates of the causal effect of feature-selection placement.

_Full repeated-CV run across 30 seeds; estimates seed/fold robustness of the leakage gap and feature stability (within-cohort, diagnostic; not biomarker discovery)._
