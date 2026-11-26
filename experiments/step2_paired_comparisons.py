"""
experiments/step2_paired_comparisons.py

2. PAIRED COMPARISONS
- Computes paired bootstrap differences (1,000 resamples of queries, seed 42) with 95% CIs.
- Evaluated on C1 for:
  - P1-dn
  - P2a-dn (0.4)
  - P2b-dn (0.4)
  - P3-dn
- Comparisons:
  - M6 - M5
  - M7 - M6
  - M7 - M5
  - M9 - M7
  - M7 - M2
  - M7 - M8
  - M6 - M2
- Metrics: MRR, hit@10
- Saves to experiments/results/paired_differences.csv.
"""

import csv
import pickle
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

RESULTS_DIR = ROOT / "experiments" / "results"

def paired_bootstrap(diff_array: np.ndarray, n_boot: int = 1000, seed: int = 42):
    n = len(diff_array)
    rng = np.random.RandomState(seed)
    boot_indices = rng.randint(0, n, size=(n_boot, n))
    
    boot_means = np.empty(n_boot, dtype=np.float64)
    for b in range(n_boot):
        boot_means[b] = np.mean(diff_array[boot_indices[b]])
        
    mean_val = float(np.mean(diff_array))
    ci_low = float(np.percentile(boot_means, 2.5))
    ci_high = float(np.percentile(boot_means, 97.5))
    return mean_val, ci_low, ci_high

def main():
    cache_path = RESULTS_DIR / "query_metrics_dn.pkl"
    if not cache_path.exists():
        print(f"Error: Cache file {cache_path} not found. Run step1_duplicate_neutral_run.py first.")
        sys.exit(1)

    with open(cache_path, "rb") as f:
        data = pickle.load(f)

    all_query_metrics = data["all_query_metrics"]

    protocols = [
        ("P1-dn", "P1_dn"),
        ("P2a-dn (0.4)", "P2a_dn_0.4"),
        ("P2b-dn (0.4)", "P2b_dn_0.4"),
        ("P3-dn", "P3_dn"),
    ]

    comparisons = [
        ("M6 - M5", "M6", "M5"),
        ("M7 - M6", "M7", "M6"),
        ("M7 - M5", "M7", "M5"),
        ("M9 - M7", "M9", "M7"),
        ("M7 - M2", "M7", "M2"),
        ("M7 - M8", "M7", "M8"),
        ("M6 - M2", "M6", "M2"),
    ]

    metrics = [
        ("mrr", "mrr"),
        ("hit@10", "hit10"),
    ]

    rows = []

    for proto_label, proto_key in protocols:
        qm = all_query_metrics[proto_key]
        n_queries = len(qm["M1"]["mrr"])

        # Precompute the exact 1000 resamples for this protocol queries
        rng = np.random.RandomState(42)
        boot_indices = rng.randint(0, n_queries, size=(1000, n_queries))

        for comp_label, mA, mB in comparisons:
            for metric_label, metric_key in metrics:
                arrA = qm[mA][metric_key]
                arrB = qm[mB][metric_key]
                diffs = arrA - arrB

                # Paired bootstrap across the same 1000 resamples
                boot_means = np.empty(1000, dtype=np.float64)
                for b in range(1000):
                    boot_means[b] = np.mean(diffs[boot_indices[b]])

                mean_diff = float(np.mean(diffs))
                ci_low = float(np.percentile(boot_means, 2.5))
                ci_high = float(np.percentile(boot_means, 97.5))

                rows.append({
                    "protocol": proto_label,
                    "comparison": comp_label,
                    "metric": metric_label,
                    "mean_diff": mean_diff,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                })

    out_csv = RESULTS_DIR / "paired_differences.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["protocol", "comparison", "metric", "mean_diff", "ci_low", "ci_high"])
        w.writeheader()
        w.writerows(rows)

    print(f"Saved {out_csv} ({len(rows)} rows)")

    # Print summary table
    print("\n=== PAIRED DIFFERENCES SUMMARY ===")
    for r in rows:
        print(f"{r['protocol']:<14} | {r['comparison']:<8} | {r['metric']:<6} | Mean: {r['mean_diff']:+.5f} (95% CI: [{r['ci_low']:+.5f}, {r['ci_high']:+.5f}])")

if __name__ == "__main__":
    main()

