"""
experiments/step3_run_all_methods.py

Executes Step 3:
Runs Methods M1 to M9 across Configurations C1, C2, C3, C4 and Protocols:
- P1
- P2a (0.3)
- P2a (0.4)
- P2a (0.5)
- P2b (0.4)
- P3 (primary ATC-restricted candidates)
- P3-all (all candidates)

For C1, also runs and saves size-stratified evaluations:
- Small (1 to 10)
- Medium (11 to 50)
- Large (>50)

Saves CSVs to experiments/results/metrics_<config>_<protocol>.csv with columns:
method, n_queries, metric, mean, ci_low, ci_high
"""

import csv
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_atc_l3, precompute_matrices,
    evaluate_protocol, bootstrap_ci, bootstrap_enrichment,
    DB_A_PATH, DB_R_PATH, GOLD_PRIMARY_PATH, DATA_DIR, RESULTS_DIR
)

PROTOCOLS = [
    ("P1", {}),
    ("P2a_0.3", {"theta": 0.3}),
    ("P2a_0.4", {"theta": 0.4}),
    ("P2a_0.5", {"theta": 0.5}),
    ("P2b_0.4", {}),
    ("P3", {}),
    ("P3_all", {}),
]

def save_metrics_csv(filename: Path, rows: list):
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "n_queries", "metric", "mean", "ci_low", "ci_high"])
        for r in rows:
            writer.writerow([
                r["method"],
                r["n_queries"],
                r["metric"],
                f"{r['mean']:.6f}",
                f"{r['ci_low']:.6f}",
                f"{r['ci_high']:.6f}",
            ])

def summarize_subset(query_metrics: dict, subset_indices: list, seed: int = 42):
    sub_n = len(subset_indices)
    if sub_n == 0:
        return []
    rows = []
    idx_arr = np.array(subset_indices, dtype=np.int32)
    for m in ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"]:
        qm = query_metrics[m]
        h10 = qm["hit10"][idx_arr]
        h20 = qm["hit20"][idx_arr]
        p10 = qm["prec10"][idx_arr]
        mrr = qm["mrr"][idx_arr]
        pct = qm["pct_rank"][idx_arr]
        rnd = qm["rand_hit10"][idx_arr]

        m_h10, l_h10, h_h10 = bootstrap_ci(h10, n_boot=1000, seed=seed)
        rows.append({"method": m, "n_queries": sub_n, "metric": "hit@10", "mean": m_h10, "ci_low": l_h10, "ci_high": h_h10})

        m_h20, l_h20, h_h20 = bootstrap_ci(h20, n_boot=1000, seed=seed)
        rows.append({"method": m, "n_queries": sub_n, "metric": "hit@20", "mean": m_h20, "ci_low": l_h20, "ci_high": h_h20})

        m_p10, l_p10, h_p10 = bootstrap_ci(p10, n_boot=1000, seed=seed)
        rows.append({"method": m, "n_queries": sub_n, "metric": "precision@10", "mean": m_p10, "ci_low": l_p10, "ci_high": h_p10})

        m_mrr, l_mrr, h_mrr = bootstrap_ci(mrr, n_boot=1000, seed=seed)
        rows.append({"method": m, "n_queries": sub_n, "metric": "mrr", "mean": m_mrr, "ci_low": l_mrr, "ci_high": h_mrr})

        m_pct, l_pct, h_pct = bootstrap_ci(pct, n_boot=1000, seed=seed)
        rows.append({"method": m, "n_queries": sub_n, "metric": "percentile_rank", "mean": m_pct, "ci_low": l_pct, "ci_high": h_pct})

        m_enr, l_enr, h_enr = bootstrap_enrichment(h10, rnd, n_boot=1000, seed=seed)
        rows.append({"method": m, "n_queries": sub_n, "metric": "enrichment", "mean": m_enr, "ci_low": l_enr, "ci_high": h_enr})
    return rows

def run_configuration(config_name: str, db_path: Path, gold_path: Path, exclude_inorganics: bool):
    print(f"\n=======================================================")
    print(f"RUNNING CONFIGURATION: {config_name}")
    print(f"Database: {db_path.name}, Gold: {gold_path.name}, Exclude inorganics: {exclude_inorganics}")
    print(f"=======================================================")

    t0 = time.time()
    drugs = load_database(db_path, exclude_inorganics=exclude_inorganics)
    pool_cids = sorted(list(drugs.keys()))
    print(f"Loaded {len(pool_cids)} pool drugs ({time.time() - t0:.2f}s)")

    drug_inds, _ = load_gold_pairs(gold_path, set(pool_cids))
    print(f"Loaded {len(drug_inds)} mapped drugs from gold")

    cid_atc = load_atc_l3(set(pool_cids))
    print(f"Loaded {sum(1 for c in pool_cids if len(cid_atc[c]) > 0)} drugs with ATC")

    t_mat = time.time()
    matrices = precompute_matrices(drugs, pool_cids)
    print(f"Precomputed S, T, E matrices ({time.time() - t_mat:.2f}s)")

    for proto_name, kwargs in PROTOCOLS:
        t_p = time.time()
        print(f"  Evaluating {config_name} - {proto_name} ...", end="", flush=True)
        summary, queries, query_metrics = evaluate_protocol(
            proto_name, pool_cids, matrices, drug_inds, cid_atc, seed=42, **kwargs
        )
        out_csv = RESULTS_DIR / f"metrics_{config_name}_{proto_name}.csv"
        save_metrics_csv(out_csv, summary)
        print(f" Done ({len(queries)} queries, {time.time() - t_p:.2f}s) -> {out_csv.name}")

        # Size strata for C1
        if config_name == "C1":
            sizes = [q["R"] for q in queries]
            small_idx = [i for i, s in enumerate(sizes) if 1 <= s <= 10]
            med_idx = [i for i, s in enumerate(sizes) if 11 <= s <= 50]
            large_idx = [i for i, s in enumerate(sizes) if s > 50]

            for s_name, s_idx in [("small", small_idx), ("medium", med_idx), ("large", large_idx)]:
                sub_rows = summarize_subset(query_metrics, s_idx, seed=42)
                sub_csv = RESULTS_DIR / f"metrics_C1_{proto_name}_{s_name}.csv"
                save_metrics_csv(sub_csv, sub_rows)

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    start_total = time.time()

    # 1. C1: DB-A Primary
    run_configuration("C1", DB_A_PATH, GOLD_PRIMARY_PATH, exclude_inorganics=True)

    # 2. C2: DB-R Primary
    run_configuration("C2", DB_R_PATH, GOLD_PRIMARY_PATH, exclude_inorganics=True)

    # 3. C3: DB-A Rule (a) only
    c3_gold = DATA_DIR / "gold_pairs_rule_a_only.csv"
    run_configuration("C3", DB_A_PATH, c3_gold, exclude_inorganics=True)

    # 4. C4: DB-A with Inorganics
    c4_gold = DATA_DIR / "gold_pairs_with_inorganics.csv"
    run_configuration("C4", DB_A_PATH, c4_gold, exclude_inorganics=False)

    print(f"\nAll configurations and protocols completed in {time.time() - start_total:.2f}s.")

if __name__ == "__main__":
    main()

