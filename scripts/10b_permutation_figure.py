#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "permutation" / "permutation_b1000_fixed_null_distributions.csv"
OUT = ROOT / "figures"


def main():
    data = pd.read_csv(DATA)
    observed = data.loc[data["permutation_id"] == 0].iloc[0]
    null = data.loc[data["permutation_id"] > 0]

    values = np.concatenate([
        null["leaky_auroc"].to_numpy(),
        null["guarded_auroc"].to_numpy(),
        np.array([0.5, observed["leaky_auroc"], observed["guarded_auroc"]]),
    ])
    bins = np.linspace(values.min() - 0.02, values.max() + 0.02, 26)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    axes[0].hist(null["leaky_auroc"], bins=bins, color="#E99695", edgecolor="#444444", alpha=0.88, label="naive null")
    axes[0].hist(null["guarded_auroc"], bins=bins, color="#C2E0C6", edgecolor="#444444", alpha=0.88, label="guarded null")
    axes[0].axvline(0.5, color="#555555", linestyle=":", linewidth=1.5, label="chance")
    axes[0].axvline(observed["leaky_auroc"], color="#C53030", linewidth=2, label="observed naive")
    axes[0].axvline(observed["guarded_auroc"], color="#1F7A3D", linewidth=2, label="observed guarded")
    axes[0].set(title="Fixed-direction AUROC under shuffled labels", xlabel="AUROC", ylabel="Frequency")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].hist(null["gap_auroc"], bins=15, color="#D9E2F3", edgecolor="#444444")
    axes[1].axvline(0, color="#555555", linestyle=":", linewidth=1.5)
    axes[1].axvline(observed["gap_auroc"], color="#2B6CB0", linewidth=2, label="observed contrast")
    axes[1].set(title="Null workflow contrast (naive − guarded)", xlabel="ΔAUROC", ylabel="Frequency")
    axes[1].legend(frameon=False, fontsize=8)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.savefig(OUT / "permutation_b1000_fixed_null.png", dpi=300)
    fig.savefig(OUT / "permutation_b1000_fixed_null.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
