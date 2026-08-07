# Evidence Audit Dashboard — Notes

The six-panel dashboard is generated from committed outputs only by `scripts/11_plot_evidence_audit_dashboard.R`.

1. Whole-workflow contrast: ΔAUROC 0.0799 (95% CI 0.0103–0.1552; p=0.028) and ΔPR-AUC 0.0789 (95% CI −0.0426–0.1778; p=0.191).
2. Fixed-orientation B=1000 permutation control: naive null AUROC mean 0.8778; guarded null AUROC mean 0.4917.
3. Thirty-seed descriptive contrast: AUROC is positive in 28/30 seeds; PR-AUC in 22/30.
4. Feature stability: submission-aligned single-run Nogueira 0.5128 and mean Jaccard 0.3487, with the 30-seed distributions shown.
5. Same-study-family transportability: guarded discovery AUROC 0.7032 versus GSE25065 AUROC 0.5633.
6. Reproducibility status and interpretation guardrails.

The dashboard labels the naive-versus-guarded comparison as a descriptive whole-workflow contrast because the two arms also differ in tuning and resampling. It is a methodology/audit visualization, not a biomarker claim.

Regenerate with `Rscript scripts/11_plot_evidence_audit_dashboard.R`.
