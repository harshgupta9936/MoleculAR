"""
experiments/step5_failure_analysis.py

5. FAILURE ANALYSIS (data for discussion)
For M7 on P2b-dn (0.4) and P3-dn:
- Lists 20 queries with best first-relevant rank and 20 with worst.
- For each:
  - query_cid, query_name, n_relevant
  - first_rel_cid, first_rel_name, first_rel_rank
  - S, T, E values
  - query_has_targets, query_has_se, cand_has_targets, cand_has_se
- Saves as experiments/results/failure_analysis_p2b_dn.csv and failure_analysis_p3_dn.csv.
"""

import csv
import pickle
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, precompute_matrices, compute_average_ranks,
    DB_A_PATH, RESULTS_DIR
)

def analyze_protocol(proto_label, proto_key, queries, drugs, matrices):
    print(f"\n=== FAILURE ANALYSIS: {proto_label} ===")
    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    penalties = matrices["penalties"]
    cid_to_idx = matrices["cid_to_idx"]

    has_target = {cid: len(drugs[cid]["targets"]) > 0 for cid in drugs}
    has_se = {cid: len(drugs[cid]["side_effects"]) > 0 for cid in drugs}

    query_details = []

    for q in queries:
        q_cid = q["q_cid"]
        q_idx = q["q_idx"]
        cand_idx = q["cand_indices"]
        rel_idx = q["rel_indices"]

        s_row = S[q_idx]
        t_row = T[q_idx]
        e_row = E[q_idx]

        scores = (0.4 * s_row + 0.4 * t_row + 0.2 * e_row) / penalties
        rank_map = compute_average_ranks(scores, cand_idx)

        # Find first relevant drug (min rank)
        best_cand_idx = None
        min_rank = 1e9
        best_score = -1.0

        for c_i in rel_idx:
            r = rank_map[c_i]
            sc = scores[c_i]
            if r < min_rank or (r == min_rank and sc > best_score):
                min_rank = r
                best_cand_idx = c_i
                best_score = sc

        first_rel_cid = list(cid_to_idx.keys())[list(cid_to_idx.values()).index(best_cand_idx)]

        query_details.append({
            "query_cid": q_cid,
            "query_name": drugs[q_cid]["name"],
            "n_relevant": len(rel_idx),
            "first_rel_cid": first_rel_cid,
            "first_rel_name": drugs[first_rel_cid]["name"],
            "first_rel_rank": float(min_rank),
            "S": float(S[q_idx, best_cand_idx]),
            "T": float(T[q_idx, best_cand_idx]),
            "E": float(E[q_idx, best_cand_idx]),
            "query_has_targets": bool(has_target[q_cid]),
            "query_has_se": bool(has_se[q_cid]),
            "cand_has_targets": bool(has_target[first_rel_cid]),
            "cand_has_se": bool(has_se[first_rel_cid]),
        })

    # Sort by first_rel_rank ascending (for ties in best, sort by n_relevant ascending, or stable)
    # Best 20: smallest rank
    # Worst 20: largest rank
    sorted_queries = sorted(query_details, key=lambda x: (x["first_rel_rank"], -x["n_relevant"]))
    best_20 = sorted_queries[:20]
    for r in best_20:
        r["category"] = "best_20"

    worst_sorted = sorted(query_details, key=lambda x: (-x["first_rel_rank"], x["n_relevant"]))
    worst_20 = worst_sorted[:20]
    for r in worst_20:
        r["category"] = "worst_20"

    combined = best_20 + worst_20

    out_csv = RESULTS_DIR / f"failure_analysis_{proto_key}.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "category", "query_cid", "query_name", "n_relevant",
            "first_rel_cid", "first_rel_name", "first_rel_rank",
            "S", "T", "E",
            "query_has_targets", "query_has_se", "cand_has_targets", "cand_has_se"
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(combined)

    print(f"Saved {out_csv} (40 rows)")

    print("\nTop 5 Best Queries:")
    for r in best_20[:5]:
        print(f"  {r['query_name']} -> {r['first_rel_name']} | Rank: {r['first_rel_rank']:.1f} | S={r['S']:.3f}, T={r['T']:.3f}, E={r['E']:.3f}")

    print("\nTop 5 Worst Queries:")
    for r in worst_20[:5]:
        print(f"  {r['query_name']} -> {r['first_rel_name']} | Rank: {r['first_rel_rank']:.1f} | S={r['S']:.3f}, T={r['T']:.3f}, E={r['E']:.3f}")

def main():
    cache_path = RESULTS_DIR / "query_metrics_dn.pkl"
    if not cache_path.exists():
        print(f"Error: Cache file {cache_path} not found. Run step1_duplicate_neutral_run.py first.")
        sys.exit(1)

    with open(cache_path, "rb") as f:
        data = pickle.load(f)

    all_queries = data["all_queries"]
    pool_cids = data["pool_cids"]

    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    matrices = precompute_matrices(drugs, pool_cids)

    analyze_protocol("P2b-dn (0.4)", "p2b_dn", all_queries["P2b_dn_0.4"], drugs, matrices)
    analyze_protocol("P3-dn", "p3_dn", all_queries["P3_dn"], drugs, matrices)

if __name__ == "__main__":
    main()
