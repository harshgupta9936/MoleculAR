"""
experiments/step1_duplicate_neutral_run.py

1. DUPLICATE AUDIT AND DUPLICATE-NEUTRAL RE-RUN
- Audits same-drug duplicates in the primary pool (N=2,192).
- Saves experiments/results/duplicate_pairs.csv and experiments/logs/duplicate_audit.log.
- Evaluates Methods M1 to M9 across:
  - P1-dn
  - P2a-dn (0.4)
  - P2b-dn (0.4)
  - P3-dn
- Saves metrics CSVs in experiments/results/.
- Saves detailed query metrics in experiments/results/query_metrics_dn.pkl for downstream tasks.
"""

import csv
import pickle
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_atc_l3, precompute_matrices,
    load_duplicates, evaluate_protocol,
    DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

LOGS_DIR = ROOT / "experiments" / "logs"

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== 1. DUPLICATE AUDIT AND DUPLICATE-NEUTRAL EVALUATION ===")
    t0 = time.time()

    # 1. Load Primary Database
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    cid_atc_l3 = load_atc_l3(set(pool_cids))

    print(f"Loaded primary pool: N = {len(pool_cids)} drugs")

    # 2. Duplicate Audit
    duplicates_of, pair_reasons = load_duplicates(pool_cids, DB_A_PATH)
    drugs_with_dups = [c for c, dups in duplicates_of.items() if len(dups) > 0]
    print(f"Pool drugs with >= 1 duplicate: {len(drugs_with_dups)}")
    print(f"Total duplicate pairs: {len(pair_reasons)}")

    # Save duplicate pairs CSV
    dup_pairs_path = RESULTS_DIR / "duplicate_pairs.csv"
    with open(dup_pairs_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["chembl_id_1", "name_1", "chembl_id_2", "name_2", "reasons"])
        for (c1, c2), reasons in sorted(pair_reasons.items()):
            w.writerow([c1, drugs[c1]["name"], c2, drugs[c2]["name"], "; ".join(reasons)])
    print(f"Saved {dup_pairs_path}")

    # Check query counts affected
    # In P1:
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}

    # Query counts for P1, P2a, P2b, P3
    # Check how many queries have duplicates removed
    affected_queries = {}
    total_queries_orig = {}
    total_queries_dn = {}

    matrices = precompute_matrices(drugs, pool_cids)

    protocols = [
        ("P1", "P1_dn"),
        ("P2a_0.4", "P2a_dn_0.4"),
        ("P2b_0.4", "P2b_dn_0.4"),
        ("P3", "P3_dn"),
    ]

    all_query_metrics = {}
    all_queries = {}

    for orig_name, dn_name in protocols:
        print(f"\n--- Running Duplicate-Neutral Protocol: {dn_name} ---")
        summary_rows, queries, qm = evaluate_protocol(
            protocol_name=dn_name,
            pool_cids=pool_cids,
            matrices=matrices,
            drug_inds=drug_inds,
            cid_atc_l3=cid_atc_l3,
            seed=42,
            duplicates_of=duplicates_of,
        )

        all_query_metrics[dn_name] = qm
        all_queries[dn_name] = queries

        out_csv = RESULTS_DIR / f"metrics_C1_{dn_name}.csv"
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["method", "n_queries", "metric", "mean", "ci_low", "ci_high"])
            w.writeheader()
            w.writerows(summary_rows)
        print(f"Saved {out_csv} ({len(summary_rows)} rows, {len(queries)} queries)")

    # Save detailed query metrics for downstream paired comparisons and ablation
    cache_path = RESULTS_DIR / "query_metrics_dn.pkl"
    with open(cache_path, "wb") as f:
        pickle.dump({
            "all_query_metrics": all_query_metrics,
            "all_queries": all_queries,
            "pool_cids": pool_cids,
            "duplicates_of": duplicates_of,
        }, f)
    print(f"Saved cached query metrics to {cache_path}")

    # Write audit log
    audit_log_path = LOGS_DIR / "duplicate_audit.log"
    with open(audit_log_path, "w", encoding="utf-8") as f:
        f.write("=== SAME-DRUG DUPLICATE AUDIT ===\n")
        f.write(f"Primary pool size: {len(pool_cids)}\n")
        f.write(f"Number of pool drugs with >= 1 duplicate: {len(drugs_with_dups)}\n")
        f.write(f"Total duplicate pairs: {len(pair_reasons)}\n")
        
        # Check P1 queries affected
        p1_q_with_dup_cand = [q for q in all_queries["P1_dn"] if len(duplicates_of.get(q["q_cid"], set())) > 0]
        f.write(f"P1-dn queries evaluated: {len(all_queries['P1_dn'])}\n")
        f.write(f"P1-dn queries with >= 1 duplicate in pool: {len(p1_q_with_dup_cand)}\n")
        f.write(f"P2a-dn (0.4) queries evaluated: {len(all_queries['P2a_dn_0.4'])}\n")
        f.write(f"P2b-dn (0.4) queries evaluated: {len(all_queries['P2b_dn_0.4'])}\n")
        f.write(f"P3-dn queries evaluated: {len(all_queries['P3_dn'])}\n\n")

        f.write("--- 20 Example Duplicate Pairs ---\n")
        for idx, ((c1, c2), reasons) in enumerate(sorted(pair_reasons.items())[:20]):
            f.write(f"{idx+1}. {c1} ({drugs[c1]['name']}) <-> {c2} ({drugs[c2]['name']}) | Reasons: {', '.join(reasons)}\n")

    print(f"Saved audit log to {audit_log_path}")
    print(f"Total time: {time.time() - t0:.2f}s")

if __name__ == "__main__":
    main()

