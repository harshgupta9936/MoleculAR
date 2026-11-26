"""
experiments/step4_hub_ablation.py

4. HUB PENALTY ABLATION
- Evaluates score = (0.4*S + 0.4*T + 0.2*E) / (1 + (p/100)^k) for k in {0, 1, 2, 3, 4}.
- Protocols:
  - P1-dn
  - P2a-dn (0.4)
  - P2b-dn (0.4)
  - P3-dn
- Metrics:
  - MRR, hit@10, and enrichment with 95% bootstrap CIs (1,000 resamples, seed 42).
- Recurrence metrics:
  - Mean hubness percentile of top-10 results across all queries.
  - Maximum share of top-10 appearances by any candidate.
  - Top-1% candidates' share of total top-10 appearances.
- Saves to experiments/results/hub_penalty_ablation.csv.
"""

import csv
import pickle
import sys
import time
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, precompute_matrices, compute_average_ranks,
    bootstrap_ci, bootstrap_enrichment, DB_A_PATH, RESULTS_DIR
)

def evaluate_hub_ablation():
    cache_path = RESULTS_DIR / "query_metrics_dn.pkl"
    if not cache_path.exists():
        print(f"Error: Cache file {cache_path} not found. Run step1_duplicate_neutral_run.py first.")
        sys.exit(1)

    with open(cache_path, "rb") as f:
        data = pickle.load(f)

    all_queries = data["all_queries"]
    pool_cids = data["pool_cids"]
    N = len(pool_cids)

    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    matrices = precompute_matrices(drugs, pool_cids)

    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    hub_pcts = matrices["hub_pcts"]

    # Raw composite score (0.4S + 0.4T + 0.2E)
    raw_score_matrix = 0.4 * S + 0.4 * T + 0.2 * E

    protocols = [
        ("P1-dn", "P1_dn"),
        ("P2a-dn (0.4)", "P2a_dn_0.4"),
        ("P2b-dn (0.4)", "P2b_dn_0.4"),
        ("P3-dn", "P3_dn"),
    ]

    k_values = [0, 1, 2, 3, 4]

    results_table = []

    print("=== 4. HUB PENALTY ABLATION ===")
    t0 = time.time()

    for proto_label, proto_key in protocols:
        queries = all_queries[proto_key]
        n_queries = len(queries)
        print(f"\n--- Protocol: {proto_label} ({n_queries} queries) ---")

        for k in k_values:
            # Denominator: 1 + (p/100)^k
            if k == 0:
                penalties = np.full(N, 2.0, dtype=np.float64)
            else:
                penalties = 1.0 + (hub_pcts / 100.0) ** k

            query_mrrs = np.empty(n_queries, dtype=np.float64)
            query_hit10s = np.empty(n_queries, dtype=np.float64)
            query_rand10s = np.empty(n_queries, dtype=np.float64)

            top10_hubness_list = []
            candidate_top10_counter = Counter()

            for q_num, q in enumerate(queries):
                q_idx = q["q_idx"]
                cand_idx = q["cand_indices"]
                rel_idx = q["rel_indices"]
                M = q["M"]
                R = q["R"]

                # Random expectation for this query
                from experiments.eval_engine import rand_hit_at_k
                rand10 = rand_hit_at_k(M, R, 10)
                query_rand10s[q_num] = rand10

                scores = raw_score_matrix[q_idx] / penalties

                # Average ranks for relevance metrics
                rank_map = compute_average_ranks(scores, cand_idx)
                rel_ranks = [rank_map[idx] for idx in rel_idx]
                r_min = np.min(rel_ranks)

                query_mrrs[q_num] = 1.0 / r_min
                query_hit10s[q_num] = 1.0 if r_min <= 10.0 else 0.0

                # Top-10 candidates for recurrence metrics
                # Sort candidates descending by score
                cand_scores = scores[cand_idx]
                order = np.argsort(-cand_scores)
                top10_cand_indices = cand_idx[order[:10]]

                for c_idx in top10_cand_indices:
                    top10_hubness_list.append(hub_pcts[c_idx])
                    candidate_top10_counter[c_idx] += 1

            # Compute bootstrap CIs (seed 42)
            m_mrr, low_mrr, high_mrr = bootstrap_ci(query_mrrs, n_boot=1000, seed=42)
            m_hit10, low_hit10, high_hit10 = bootstrap_ci(query_hit10s, n_boot=1000, seed=42)
            m_enr, low_enr, high_enr = bootstrap_enrichment(query_hit10s, query_rand10s, n_boot=1000, seed=42)

            # Recurrence metrics
            mean_top10_hubness = float(np.mean(top10_hubness_list))
            total_top10_slots = len(top10_hubness_list) # 10 * n_queries

            # Max share by any single candidate
            max_cand_count = candidate_top10_counter.most_common(1)[0][1] if candidate_top10_counter else 0
            max_share = max_cand_count / total_top10_slots

            # Top-1% candidates' share of top-10 slots
            # 1% of candidate pool: ceil(0.01 * N) = 22 candidates
            n_top_1pct = int(np.ceil(0.01 * N))
            top_1pct_counts = sum(cnt for _, cnt in candidate_top10_counter.most_common(n_top_1pct))
            top1pct_share = top_1pct_counts / total_top10_slots

            row = {
                "protocol": proto_label,
                "k": k,
                "n_queries": n_queries,
                "mrr": m_mrr,
                "mrr_ci_low": low_mrr,
                "mrr_ci_high": high_mrr,
                "hit10": m_hit10,
                "hit10_ci_low": low_hit10,
                "hit10_ci_high": high_hit10,
                "enrichment": m_enr,
                "enrichment_ci_low": low_enr,
                "enrichment_ci_high": high_enr,
                "mean_top10_hubness": mean_top10_hubness,
                "max_cand_count": max_cand_count,
                "max_share": max_share,
                "top1pct_share": top1pct_share,
            }
            results_table.append(row)

            print(f"k={k}: MRR={m_mrr:.4f} [{low_mrr:.4f}, {high_mrr:.4f}] | "
                  f"Hit@10={m_hit10:.4f} [{low_hit10:.4f}, {high_hit10:.4f}] | "
                  f"Enrich={m_enr:.2f} | MeanHub={mean_top10_hubness:.1f}% | "
                  f"MaxShare={max_share*100:.2f}% | Top1%Share={top1pct_share*100:.2f}%")

    # Save to CSV
    out_csv = RESULTS_DIR / "hub_penalty_ablation.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "protocol", "k", "n_queries",
            "mrr", "mrr_ci_low", "mrr_ci_high",
            "hit10", "hit10_ci_low", "hit10_ci_high",
            "enrichment", "enrichment_ci_low", "enrichment_ci_high",
            "mean_top10_hubness", "max_cand_count", "max_share", "top1pct_share"
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results_table)

    print(f"\nSaved {out_csv} ({len(results_table)} rows)")
    print(f"Total time: {time.time() - t0:.1f}s")

if __name__ == "__main__":
    evaluate_hub_ablation()
