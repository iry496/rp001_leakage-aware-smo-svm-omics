# Selector K-sweep (full run, GSE25055 only)

- Selector fixed (Welch t-test top-K); K in {25, 50, 100, 200}. Seed 20260620.
- Leaky: global FS top-K before 5x5 repeated CV (cost 1). Guarded: nested 5-outer x 5-inner; FS + cost tuning in training folds.
- K=100 is the anchor and reproduces the submission-aligned run (naive ~0.7830, guarded ~0.7032, Nogueira ~0.5128).

## Results by K
- K=25: AUROC leaky 0.8344 / guarded 0.7342 (gap 0.1001); PR-AUC leaky 0.4961 / guarded 0.3808 (gap 0.1154); bal.acc 0.5780 / 0.5346; MCC 0.2662 / 0.1246; sens 0.1825 / 0.1053; spec 0.9735 / 0.9639; Nogueira 0.3513, mean Jaccard 0.2180, stable core 3, unstable tail 44.
- K=50: AUROC leaky 0.8161 / guarded 0.7204 (gap 0.0956); PR-AUC leaky 0.4695 / guarded 0.3874 (gap 0.0821); bal.acc 0.5796 / 0.5581; MCC 0.2564 / 0.2166; sens 0.1930 / 0.1404; spec 0.9663 / 0.9759; Nogueira 0.3846, mean Jaccard 0.2413, stable core 7, unstable tail 88.
- K=100: AUROC leaky 0.7830 / guarded 0.7032 (gap 0.0799); PR-AUC leaky 0.4257 / guarded 0.3476 (gap 0.0781); bal.acc 0.5525 / 0.5724; MCC 0.1945 / 0.2138; sens 0.1298 / 0.1930; spec 0.9751 / 0.9518; Nogueira 0.5128, mean Jaccard 0.3487, stable core 26, unstable tail 105.
- K=200: AUROC leaky 0.8212 / guarded 0.7054 (gap 0.1157); PR-AUC leaky 0.4759 / guarded 0.3463 (gap 0.1296); bal.acc 0.5922 / 0.5589; MCC 0.2760 / 0.1899; sens 0.2246 / 0.1579; spec 0.9598 / 0.9598; Nogueira 0.5752, mean Jaccard 0.4087, stable core 66, unstable tail 173.

- Total runtime: 0.7 min.

_Sensitivity analysis only; within-cohort, diagnostic; not biomarker discovery._
