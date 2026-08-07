#!/usr/bin/env python3
"""scripts/10_manuscript_figures.py

Generate three conceptual/diagnostic manuscript figures for the leakage-aware
SMO/SVM audit. NO new analysis is performed and NO result values are computed
here: the performance figure reads the committed bootstrap CIs, and the
feature-frequency figure reads the committed per-fold selected-feature list.
The workflow schematic is purely conceptual (no numbers).

Outputs (PNG, 200 dpi) -> figures/manuscript/:
  - fig_workflow_schematic.png      (leaky vs guarded pipeline diagram)
  - fig_performance_forest.png      (AUROC / PR-AUC, leaky vs guarded, 95% CI)
  - fig_feature_frequency.png       (selection frequency across 5 outer folds)

Inputs (committed):
  - tables/uncertainty/bootstrap_ci.csv
  - results/pilot_gse25055/nested_selected_features_by_fold.csv
"""
import csv, collections, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = "figures/manuscript"
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- Figure 1
def workflow_schematic():
    fig, ax = plt.subplots(figsize=(11.6, 6.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    def box(x, y, w, h, text, fc, ec="#333333", fs=8.5):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.12",
                                    fc=fc, ec=ec, lw=1.2))
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, wrap=True)

    def arrow(x1, y1, x2, y2, color="#444444"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=14, lw=1.4, color=color))

    ax.text(2.5, 9.6, "Naive (global-selection) workflow", ha="center", fontsize=12, fontweight="bold", color="#B03030")
    ax.text(7.5, 9.6, "Guarded (nested) workflow", ha="center", fontsize=12, fontweight="bold", color="#1E7A3C")

    # Leaky column
    box(0.6, 8.3, 3.8, 0.9, "Full dataset (all samples + labels)", "#F2F2F2")
    box(0.6, 6.7, 3.8, 1.0, "Feature selection (t-test top-100)\non the FULL dataset", "#F4C7C3")
    box(0.6, 5.1, 3.8, 0.9, "5x5 cross-validation\n(train / test split)", "#F2F2F2")
    box(0.6, 3.7, 3.8, 0.8, "Reported performance", "#F2F2F2")
    arrow(2.5, 8.3, 2.5, 7.7); arrow(2.5, 6.7, 2.5, 6.0); arrow(2.5, 5.1, 2.5, 4.5)
    ax.text(4.55, 7.2, "leakage:\nFS sees\ntest-fold\nlabels", ha="left", va="center",
            fontsize=8, color="#B03030", fontweight="bold")

    # Guarded column
    box(5.6, 8.3, 3.8, 0.9, "Full dataset (all samples + labels)", "#F2F2F2")
    box(5.6, 6.9, 3.8, 0.9, "Outer 5-fold split\n(outer test held out)", "#F2F2F2")
    box(5.6, 5.2, 3.8, 1.1, "Inner CV on OUTER-TRAIN only:\nfeature selection + cost tuning\n(test fold never seen)", "#C6E5CE")
    box(5.6, 3.6, 3.8, 0.9, "Evaluate on untouched\nouter-test fold", "#F2F2F2")
    arrow(7.5, 8.3, 7.5, 7.8); arrow(7.5, 6.9, 7.5, 6.3); arrow(7.5, 5.2, 7.5, 4.5)
    ax.text(9.55, 5.75, "FS isolated\nfrom test\nfold", ha="left", va="center",
            fontsize=8, color="#1E7A3C", fontweight="bold")

    ax.text(5.0, 2.6, "Both workflows use the same linear SVM and top-K = 100 selector. They also differ in cost policy\n"
                          "and resampling, so this is a descriptive comparison rather than a matched one-factor ablation.",
            ha="center", va="center", fontsize=8.5, style="italic")
    plt.tight_layout()
    fig.savefig(f"{OUT}/fig_workflow_schematic.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------- Figure 2
def performance_forest():
    rows = {}
    with open("tables/uncertainty/bootstrap_ci.csv") as f:
        for r in csv.DictReader(f):
            if r["cohort"] == "GSE25055" and r["metric"] in ("auroc", "pr_auc"):
                rows[(r["pipeline"], r["metric"])] = (float(r["point"]), float(r["ci_lo"]), float(r["ci_hi"]))
    # order top->bottom
    items = [
        ("AUROC — Leaky",   rows[("A_leaky_baseline","auroc")],  "#C0392B"),
        ("AUROC — Guarded", rows[("B_guarded_nested","auroc")],  "#1E8449"),
        ("PR-AUC — Leaky",  rows[("A_leaky_baseline","pr_auc")], "#C0392B"),
        ("PR-AUC — Guarded",rows[("B_guarded_nested","pr_auc")], "#1E8449"),
    ]
    fig, ax = plt.subplots(figsize=(8.2, 3.8))
    ys = list(range(len(items)))[::-1]
    for y,(label,(p,lo,hi),c) in zip(ys, items):
        ax.errorbar(p, y, xerr=[[p-lo],[hi-p]], fmt="o", color=c, ecolor=c,
                    elinewidth=2, capsize=5, markersize=8)
        ax.text(hi+0.012, y, f"{p:.3f} [{lo:.3f}, {hi:.3f}]", va="center", fontsize=8.5)
    ax.set_yticks(ys); ax.set_yticklabels([it[0] for it in items], fontsize=9.5)
    ax.set_xlim(0.25, 1.02); ax.set_xlabel("Score (95% bootstrap CI)", fontsize=10)
    ax.axhline(1.5, color="#cccccc", lw=0.8, ls=":")
    ax.set_title("Leaky vs guarded performance on GSE25055 (95% bootstrap CIs)", fontsize=11)
    ax.grid(axis="x", color="#eeeeee")
    plt.tight_layout()
    fig.savefig(f"{OUT}/fig_performance_forest.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------- Figure 6
def feature_frequency():
    byfeat = collections.defaultdict(set)
    with open("results/pilot_gse25055/nested_selected_features_by_fold.csv") as f:
        for r in csv.DictReader(f):
            byfeat[r["feature"]].add(r["fold"])
    freq = collections.Counter(len(v) for v in byfeat.values())
    xs = [1,2,3,4,5]; ys = [freq.get(k,0) for k in xs]
    colors = ["#C0392B","#E59866","#BDC3C7","#A9CCE3","#1E8449"]
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    bars = ax.bar([str(x) for x in xs], ys, color=colors, edgecolor="#333333", lw=0.8)
    for b,v in zip(bars,ys):
        ax.text(b.get_x()+b.get_width()/2, v+1.5, str(v), ha="center", fontsize=9.5, fontweight="bold")
    ax.annotate(f"unstable tail\n({ys[0]} probes in 1 fold)", xy=(0, ys[0]), xytext=(0.6,86),
                fontsize=8.5, color="#B03030",
                arrowprops=dict(arrowstyle="->", color="#B03030"))
    ax.annotate(f"stable core\n({ys[4]} probes in all 5 folds)", xy=(4, ys[4]), xytext=(2.7,55),
                fontsize=8.5, color="#1E7A3C",
                arrowprops=dict(arrowstyle="->", color="#1E7A3C"))
    ax.set_xlabel("Number of outer folds selecting a feature (of 5)", fontsize=10)
    ax.set_ylabel("Number of features", fontsize=10)
    ax.set_title("Feature-selection frequency across outer folds (GSE25055 guarded pipeline)", fontsize=11)
    ax.set_ylim(0, max(ys)*1.2)
    ax.spines[["top","right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{OUT}/fig_feature_frequency.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return dict(freq)

if __name__ == "__main__":
    workflow_schematic()
    performance_forest()
    dist = feature_frequency()
    print("fold-frequency distribution used:", {k: dist.get(k,0) for k in (1,2,3,4,5)})
    print("figures written to", OUT)
