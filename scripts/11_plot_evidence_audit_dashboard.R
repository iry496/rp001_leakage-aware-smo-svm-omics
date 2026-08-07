#!/usr/bin/env Rscript
# =============================================================================
# scripts/11_plot_evidence_audit_dashboard.R
# -----------------------------------------------------------------------------
# Six-panel Evidence Audit dashboard from committed outputs. NO NEW MODELING.
# Panels: (1) single-seed whole-workflow contrast with CI; (2) B=1000 permutation null;
# (3) 30-seed workflow-contrast distribution; (4) feature-stability distribution;
# (5) same-study-family external-validation drop; (6) reproducibility/audit status.
# Writes figures/evidence_audit_dashboard.png and .pdf.
# =============================================================================
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
dci <- rd("tables/uncertainty/delta_auroc_prauc_ci.csv")
nd  <- rd("results/permutation/permutation_b1000_fixed_null_distributions.csv")
nd  <- nd[toupper(as.character(nd$is_identity)) == "FALSE", ]
prm <- rd("results/permutation/permutation_b1000_fixed_pvalues.csv")
gap <- rd("tables/repeated_cv/leakage_gap_by_seed.csv")
sbs <- rd("tables/repeated_cv/stability_by_seed.csv")
bci <- rd("tables/uncertainty/bootstrap_ci.csv")
ext <- rd("results/external_validation_gse25065/gse25065_external_metrics.csv")
fss <- rd("tables/pilot_gse25055/feature_stability_summary.csv")
bg <- function(p, m, c) bci[bci$pipeline == p & bci$metric == m, c][1]
pv <- function(s, c) prm[prm$statistic == s, c][1]
fs <- function(m) as.numeric(fss$value[fss$metric == m][1])

RED <- "#C0392B"; GRN <- "#1E8449"; BLU <- "#2E5A88"; GREY <- "#555555"

draw <- function() {
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

  # Panel 1: single-seed whole-workflow contrast with CI
  da <- dci[dci$metric == "delta_auroc", ]; dp <- dci[dci$metric == "delta_pr_auc", ]
  plot(NA, xlim = c(-0.1, 0.22), ylim = c(-0.5, 1.5), yaxt = "n", xlab = "contrast (naive - guarded)", ylab = "",
       main = "1. Single-seed whole-workflow contrast (95% CI)")
  axis(2, at = c(0, 1), labels = c("dPR-AUC", "dAUROC"), las = 1)
  abline(v = 0, lty = 3, col = GREY)
  for (it in list(list(y = 1, r = da), list(y = 0, r = dp))) {
    r <- it$r; points(r$point, it$y, pch = 19, col = BLU, cex = 1.4)
    segments(r$ci_lo, it$y, r$ci_hi, it$y, col = BLU, lwd = 2)
    text(r$ci_hi + 0.01, it$y, sprintf("%+.3f [%+.3f, %+.3f]\nboot p=%s", r$point, r$ci_lo, r$ci_hi, r$boot_p),
         pos = 4, cex = 0.6, xpd = NA)
  }
  text(-0.095, -0.38, "AUROC interval excludes zero; PR-AUC interval includes zero",
       pos = 4, cex = 0.55, col = GREY)

  # Panel 2: permutation null behavior
  brks <- seq(min(c(nd$leaky_auroc, nd$guarded_auroc, 0.5)) - 0.02, max(nd$leaky_auroc) + 0.02, length.out = 26)
  hL <- hist(nd$leaky_auroc, breaks = brks, plot = FALSE); hG <- hist(nd$guarded_auroc, breaks = brks, plot = FALSE)
  plot(hG, col = "#C2E0C6", xlim = range(brks), ylim = c(0, max(hL$counts, hG$counts)),
       main = "2. Permutation null, B=1000\n(diagnostic of invalid naive workflow)", xlab = "AUROC under shuffled labels")
  plot(hL, col = "#E99695", add = TRUE)
  abline(v = 0.5, lty = 3, col = GREY)
  abline(v = as.numeric(pv("leaky_auroc", "observed")), col = RED, lwd = 2)
  abline(v = as.numeric(pv("guarded_auroc", "observed")), col = GRN, lwd = 2)
  legend("topleft", fill = c("#E99695", "#C2E0C6"), legend = c("naive null", "guarded null"), bty = "n", cex = 0.7)

  # Panel 3: 30-seed workflow-contrast distribution
  npos <- sum(gap$gap_auroc > 0)
  boxplot(gap$gap_auroc, col = "#D9E2F3", ylab = "dAUROC (naive - guarded)",
          main = sprintf("3. 30-seed whole-workflow contrast\n(descriptive; %d/30 > 0)", npos))
  points(jitter(rep(1, nrow(gap)), 6), gap$gap_auroc, pch = 19, col = "#34495E", cex = 0.7)
  abline(h = 0, lty = 3, col = GREY)

  # Panel 4: feature-stability distribution
  boxplot(list(Nogueira = sbs$nogueira_stability, `mean Jaccard` = sbs$mean_jaccard),
          col = "#E8DAEF", ylab = "stability index",
          main = "4. Feature stability across 30 seeds\n(moderate feature stability)")
  points(jitter(rep(1, nrow(sbs)), 6), sbs$nogueira_stability, pch = 19, col = "#6C3483", cex = 0.6)
  points(jitter(rep(2, nrow(sbs)), 6), sbs$mean_jaccard, pch = 19, col = "#1F618D", cex = 0.6)
  segments(0.7, fs("nogueira_stability"), 1.3, fs("nogueira_stability"), lty = 2, col = "#6C3483")
  segments(1.7, fs("mean_pairwise_jaccard"), 2.3, fs("mean_pairwise_jaccard"), lty = 2, col = "#1F618D")

  # Panel 5: same-study-family external validation drop
  disc <- c(as.numeric(bg("B_guarded_nested","auroc","point")),
            as.numeric(bg("B_guarded_nested","pr_auc","point")),
            as.numeric(bg("B_guarded_nested","sensitivity","point")))
  extv <- c(ext$auroc, ext$pr_auc, ext$sensitivity)
  m <- rbind(discovery = disc, external = extv)
  bp <- barplot(m, beside = TRUE, names.arg = c("AUROC","PR-AUC","Sensitivity"), ylim = c(0, 0.9),
                col = c("#1E8449", "#7FB3A6"), las = 1,
                main = "5. Same-study-family external validation\n(transportability audit)")
  legend("topright", fill = c("#1E8449", "#7FB3A6"), legend = c("discovery (guarded CV)", "external (GSE25065)"), bty = "n", cex = 0.7)
  text(bp, c(rbind(disc, extv)) + 0.03, sprintf("%.2f", c(rbind(disc, extv))), cex = 0.6)

  # Panel 6: reproducibility / audit status (text)
  plot.new(); title(main = "6. Reproducibility & audit status")
  lines <- c("[OK] Bootstrap CIs (07_bootstrap_ci.R)",
             "[OK] Permutation control B=1000 (08)",
             "[OK] Repeated nested CV, 30 seeds (09)",
             "[OK] Evidence-audit artifact (10)",
             "[OK] Fixed seeds, committed outputs",
             "[OK] Frozen external models (C=4; no external tuning)",
             "", "Framing: within-cohort, diagnostic.", "Not biomarker discovery.")
  for (i in seq_along(lines))
    text(0.02, 0.92 - (i - 1) * 0.1, lines[i], pos = 4, cex = 0.8,
         col = if (startsWith(lines[i], "[OK]")) GRN else GREY, xpd = NA)

  mtext("Evidence Audit Dashboard | GSE25055 discovery; external audits GSE25065/GSE41998; methodology, not biomarker discovery",
        outer = TRUE, cex = 0.85, font = 2)
}

dir.create("figures", showWarnings = FALSE)
png("figures/evidence_audit_dashboard.png", width = 1350, height = 820, res = 110); draw(); dev.off()
pdf("figures/evidence_audit_dashboard.pdf", width = 13.5, height = 8.2); draw(); dev.off()
cat("[dashboard] wrote figures/evidence_audit_dashboard.{png,pdf}\n")
