# Permutation Control - b1000_fixed run (GSE25055 only)

- Permutations: identity + 1000 shuffled. Seed 20260620 (label shuffles offset +100000).
- Fast Welch selector: TRUE ; parallel: TRUE (4 workers).
- Folds regenerated per permutation (stratified on shuffled labels); FS rerun per permutation.
- Leaky arm: global t-test top-100 before CV. Guarded arm: t-test inside training folds only.

## Identity reproduction check
- leaky AUROC = 0.7830 (ref 0.7830); guarded AUROC = 0.7032 (ref 0.7032). Within 0.01: PASS.

## Selector validation (original vs fast Welch)
- identity: overlap 100/100, exact_order=TRUE.
- perm1: overlap 100/100, exact_order=TRUE.
- perm2: overlap 100/100, exact_order=TRUE.
- perm3: overlap 100/100, exact_order=TRUE.
- perm4: overlap 100/100, exact_order=TRUE.
- perm5: overlap 100/100, exact_order=TRUE.

## Null summary (shuffled labels)
- leaky AUROC null: mean 0.8778 (0.6505-0.9760); frac > 0.5 = 1.000.
- guarded AUROC null: mean 0.4917 (0.3523-0.6342); frac > 0.5 = 0.430.
- leakage gap null: mean 0.3860 (0.1408-0.5967).

## Runtime
- mean compute 9.2 s/perm; wall 2.4 s/perm with 4 workers; total 39.5 min.
- old (serial, original selector) was 361.8 s/perm.
- estimated B=200: 7.9 min (0.13 h) ; B=1000: 39.5 min (0.66 h) (at current wall throughput).

_This permutation run estimates null behavior under shuffled labels. The leaky null distribution should be interpreted as a diagnostic of feature-selection leakage, not as evidence of biological signal. It supports the interpretation that feature-selection leakage can inflate apparent performance._
