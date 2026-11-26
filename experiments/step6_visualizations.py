"""
experiments/step6_visualizations.py

Generates evaluation figures:
1. figures/mrr_comparison.png: Grouped bar chart of MRR across M1 to M9 for P1, P2a (0.4), and P3 with 95% CI error bars.
2. figures/enrichment_comparison.png: Grouped bar chart of Enrichment across M1 to M9 for P1, P2a (0.4), and P3 with 95% CI error bars.
"""

import csv
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "experiments" / "results"
FIGURES_DIR = ROOT / "experiments" / "figures"

def load_metric_data(csv_path: Path, target_metric: str):
    data = {}
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            if r["metric"] == target_metric:
                m = r["method"]
                mean = float(r["mean"])
                low = float(r["ci_low"])
                high = float(r["ci_high"])
                data[m] = (mean, low, high)
    return data

def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    methods = ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"]

    protocols = [
        ("P1 (All)", RESULTS_DIR / "metrics_C1_P1.csv", "#1f77b4"),
        ("P2a (theta=0.4)", RESULTS_DIR / "metrics_C1_P2a_0.4.csv", "#ff7f0e"),
        ("P3 (Cross-Class ATC)", RESULTS_DIR / "metrics_C1_P3.csv", "#2ca02c"),
    ]

    # 1. MRR Bar Chart
    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)
    x = np.arange(len(methods))
    width = 0.25

    for idx, (p_label, csv_file, color) in enumerate(protocols):
        p_data = load_metric_data(csv_file, "mrr")
        means = [p_data[m][0] for m in methods]
        lows = [p_data[m][1] for m in methods]
        highs = [p_data[m][2] for m in methods]
        err_below = [m - l for m, l in zip(means, lows)]
        err_above = [h - m for m, h in zip(means, highs)]

        offset = (idx - 1) * width
        ax.bar(
            x + offset, means, width, yerr=[err_below, err_above],
            label=p_label, color=color, capsize=4, alpha=0.85, edgecolor="black", linewidth=0.5
        )

    ax.set_ylabel("Mean Reciprocal Rank (MRR)", fontsize=12)
    ax.set_title("MRR Comparison Across Methods (DB-A Primary Pool, 95% Bootstrap CIs)", fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()

    mrr_fig_path = FIGURES_DIR / "mrr_comparison.png"
    plt.savefig(mrr_fig_path, dpi=300)
    plt.close()
    print(f"Saved MRR figure to {mrr_fig_path}")

    # 2. Enrichment Bar Chart
    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)
    for idx, (p_label, csv_file, color) in enumerate(protocols):
        p_data = load_metric_data(csv_file, "enrichment")
        means = [p_data[m][0] for m in methods]
        lows = [p_data[m][1] for m in methods]
        highs = [p_data[m][2] for m in methods]
        err_below = [m - l for m, l in zip(means, lows)]
        err_above = [h - m for m, h in zip(means, highs)]

        offset = (idx - 1) * width
        ax.bar(
            x + offset, means, width, yerr=[err_below, err_above],
            label=p_label, color=color, capsize=4, alpha=0.85, edgecolor="black", linewidth=0.5
        )

    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1.5, label="Random Baseline (1.0x)")
    ax.set_ylabel("Enrichment (Hit@10 / Random Hit@10)", fontsize=12)
    ax.set_title("Enrichment Comparison Across Methods (DB-A Primary Pool, 95% Bootstrap CIs)", fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=11)
    ax.legend(fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()

    enrich_fig_path = FIGURES_DIR / "enrichment_comparison.png"
    plt.savefig(enrich_fig_path, dpi=300)
    plt.close()
    print(f"Saved Enrichment figure to {enrich_fig_path}")

    # 3. Weight Sensitivity Heatmap
    grid_csv_path = RESULTS_DIR / "weight_grid.csv"
    if grid_csv_path.exists():
        import pandas as pd
        df_grid = pd.read_csv(grid_csv_path)
        heatmap_matrix = np.full((11, 11), np.nan)
        for _, r in df_grid.iterrows():
            i = int(round(r["w_t"] * 10))
            j = int(round(r["w_s"] * 10))
            heatmap_matrix[i, j] = r["mrr_tune"]

        best_idx = df_grid["mrr_tune"].idxmax()
        best_row = df_grid.loc[best_idx]
        best_ws = float(best_row["w_s"])
        best_wt = float(best_row["w_t"])
        best_we = float(best_row["w_e"])
        best_mrr = float(best_row["mrr_tune"])

        plt.figure(figsize=(8, 7), dpi=300)
        cmap = plt.cm.viridis
        im = plt.imshow(heatmap_matrix, origin="lower", extent=[-0.05, 1.05, -0.05, 1.05], cmap=cmap)
        cbar = plt.colorbar(im)
        cbar.set_label("MRR on Tuning Set (P2a $\\theta=0.4$)", fontsize=11)
        plt.scatter([best_ws], [best_wt], color="red", marker="*", s=200, label=f"Best ({best_ws:.1f}, {best_wt:.1f}, {best_we:.1f}): {best_mrr:.4f}")
        plt.scatter([0.4], [0.4], color="cyan", marker="o", s=80, edgecolors="black", label="Default (0.4, 0.4, 0.2)")
        plt.title("Weight Sensitivity Grid Search on P2a ($\\theta=0.4$, Tune Set)", fontsize=12, fontweight="bold")
        plt.xlabel("Structural Weight ($w_s$)", fontsize=11)
        plt.ylabel("Target Weight ($w_t$)", fontsize=11)
        plt.xticks(np.arange(0.0, 1.1, 0.1))
        plt.yticks(np.arange(0.0, 1.1, 0.1))
        plt.grid(color="white", linestyle="--", linewidth=0.5, alpha=0.5)
        plt.legend(loc="upper right", framealpha=0.9)
        plt.tight_layout()

        heatmap_path = FIGURES_DIR / "weight_heatmap.png"
        plt.savefig(heatmap_path, dpi=300)
        plt.close()
        print(f"Saved Heatmap figure to {heatmap_path}")

if __name__ == "__main__":
    main()

