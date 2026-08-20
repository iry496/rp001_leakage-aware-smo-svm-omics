# Major Output Delivery Manifest

This package contains the final working-branch outputs for the matched
feature-selection-placement ablation, unified frozen-model external validation,
and the revised NCI ODS poster abstract.

## Included analysis code

- `R/metrics.R`
- `scripts/16_matched_ablation_gse25055.R`
- `scripts/17_external_validation_unified.R`

## Included matched-ablation outputs

- Final 30-repeat predictions, fold assignments, tuning records, selected
  features, repeat-averaged predictions, design contract, paired metrics,
  paired summaries, conditional patient-bootstrap intervals, notes, and
  publication-ready PDF/PNG figures.
- Primary result: mean AUROC 0.7870 versus 0.7256; paired mean difference
  0.0615; positive in 30/30 repeats; empirical 2.5th–97.5th percentile
  split-distribution interval 0.0084–0.1201.
- The split-distribution interval is a robustness distribution across
  correlated repeated partitions, not an independence-based confidence
  interval.

## Included external-validation outputs

- Single-fit frozen-model manifest, feature manifest, full-precision external
  predictions, point estimates, conditional bootstrap intervals, cohort
  context, regression checks, notes, and PDF/PNG figures.
- GSE25065 AUROC: 0.5633; conditional 95% bootstrap interval 0.4643–0.6684.
- GSE41998 primary AUROC: 0.6922; conditional 95% bootstrap interval
  0.6146–0.7599.
- Bootstrap intervals are conditional on the frozen model and stored external
  predictions; the discovery workflow is not refit inside bootstrap samples.

## Included abstract

- `deliverables/NCI_ODS_Abstract_266_words.md`
- The abstract contains no DOI.

Temporary smoke outputs, processed expression matrices, raw data, RDS caches,
and unrelated archived manuscript files are excluded.
