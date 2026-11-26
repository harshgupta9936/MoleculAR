"""
experiments/step3_weight_stability.py

3. WEIGHT SENSITIVITY, DONE PROPERLY
- 66 triples (w_s, w_t, w_e) with step 0.1 summing to 1.0.
- Evaluated on P2b-dn (0.4) and, separately, on P1-dn.
- 20 random 50/50 splits of query drugs (seeds 0 to 19).
- For each split:
  - Best tune triple targeting MRR.
  - Test MRR and test hit@10 for the best tune triple.
  - Test MRR and test hit@10 for default triple (0.4, 0.4, 0.2).
  - Default triple's rank among the 66 triples on tune set.
- Reports selection frequency across 20 splits, and mean/std of test difference (best minus default).
- Saves to experiments/results/weight_stability_p2b_dn.csv and weight_stability_p1_dn.csv.
"""

import csv
import pickle
import sys
import time
from pathlib import Path
from collections import Counter
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_atc_l3, precompute_matrices,
    compute_average_ranks, DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

def generate_weight_triples():
    triples = []
    for i in range(11):
        for j in range(11 - i):
            k = 10 - i - j
            w_s = round(i / 10.0, 1)
            w_t = round(j / 10.0, 1)
            w_e = round(k / 10.0, 1)
            triples.append((w_s, w_t, w_e))
    return triples

def evaluate_objective(objective_name, queries, matrices, triples, default_triple=(0.4, 0.4, 0.2)):
    print(f"\n==========================================")
    print(f"RUNNING WEIGHT STABILITY: {objective_name}")
    print(f"Total queries: {len(queries)}, Triples: {len(triples)}")
    print(f"==========================================")

    t0 = time.time()
    n_queries = len(queries)
    n_triples = len(triples)

    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    penalties = matrices["penalties"]

    # Precompute per-query MRR and Hit@10 for all 66 triples
    # query_mrr[q, t], query_hit10[q, t]
    query_mrr = np.empty((n_queries, n_triples), dtype=np.float64)
    query_hit10 = np.empty((n_queries, n_triples), dtype=np.float64)

    default_idx = triples.index(default_triple)

    print("Precomputing per-query metrics for 66 triples across all queries...")
    for q_num, q in enumerate(queries):
        q_idx = q["q_idx"]
        cand_idx = q["cand_indices"]
        rel_idx = q["rel_indices"]

        s_row = S[q_idx]
        t_row = T[q_idx]
        e_row = E[q_idx]

        for t_idx, (w_s, w_t, w_e) in enumerate(triples):
            scores = (w_s * s_row + w_t * t_row + w_e * e_row) / penalties
            # compute average ranks
            rank_map = compute_average_ranks(scores, cand_idx)
            rel_ranks = [rank_map[idx] for idx in rel_idx]
            r_min = np.min(rel_ranks)

            query_mrr[q_num, t_idx] = 1.0 / r_min
            query_hit10[q_num, t_idx] = 1.0 if r_min <= 10.0 else 0.0

        if (q_num + 1) % 400 == 0 or (q_num + 1) == n_queries:
            print(f"  Processed {q_num + 1}/{n_queries} queries ({time.time() - t0:.1f}s)")

    # Now evaluate across 20 random 50/50 splits (seeds 0 to 19)
    split_rows = []
    best_triples_list = []
    test_diffs_mrr = []
    test_diffs_hit10 = []

    for seed in range(20):
        rng = np.random.RandomState(seed)
        perm = rng.permutation(n_queries)
        half = n_queries // 2
        tune_idx = perm[:half]
        test_idx = perm[half:]

        # Mean MRR on tune set
        tune_mrrs = np.mean(query_mrr[tune_idx], axis=0) # shape (66,)
        
        # Best tune triple (highest tune MRR)
        best_t_idx = int(np.argmax(tune_mrrs))
        best_triple = triples[best_t_idx]
        best_triples_list.append(best_triple)

        # Default rank on tune set (1-based, average rank for ties)
        def_tune_mrr = tune_mrrs[default_idx]
        n_gt = np.sum(tune_mrrs > def_tune_mrr)
        n_eq = np.sum(tune_mrrs == def_tune_mrr)
        def_tune_rank = 1.0 + n_gt + 0.5 * (n_eq - 1)

        # Performance on test set
        test_mrrs = np.mean(query_mrr[test_idx], axis=0)
        test_hit10s = np.mean(query_hit10[test_idx], axis=0)

        best_test_mrr = float(test_mrrs[best_t_idx])
        best_test_hit10 = float(test_hit10s[best_t_idx])
        def_test_mrr = float(test_mrrs[default_idx])
        def_test_hit10 = float(test_hit10s[default_idx])

        diff_mrr = best_test_mrr - def_test_mrr
        diff_hit10 = best_test_hit10 - def_test_hit10

        test_diffs_mrr.append(diff_mrr)
        test_diffs_hit10.append(diff_hit10)

        split_rows.append({
            "split_seed": seed,
            "best_w_s": best_triple[0],
            "best_w_t": best_triple[1],
            "best_w_e": best_triple[2],
            "best_tune_mrr": float(tune_mrrs[best_t_idx]),
            "best_test_mrr": best_test_mrr,
            "best_test_hit10": best_test_hit10,
            "default_tune_mrr": float(def_tune_mrr),
            "default_test_mrr": def_test_mrr,
            "default_test_hit10": def_test_hit10,
            "default_tune_rank": float(def_tune_rank),
            "test_diff_mrr": float(diff_mrr),
            "test_diff_hit10": float(diff_hit10),
        })

    # Save to CSV
    out_csv = RESULTS_DIR / f"weight_stability_{objective_name}.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "split_seed", "best_w_s", "best_w_t", "best_w_e",
            "best_tune_mrr", "best_test_mrr", "best_test_hit10",
            "default_tune_mrr", "default_test_mrr", "default_test_hit10",
            "default_tune_rank", "test_diff_mrr", "test_diff_hit10",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(split_rows)
    print(f"Saved {out_csv}")

    # Compute top-line summary stats
    mean_diff_mrr = float(np.mean(test_diffs_mrr))
    std_diff_mrr = float(np.std(test_diffs_mrr, ddof=1)) if len(test_diffs_mrr) > 1 else 0.0
    mean_diff_hit10 = float(np.mean(test_diffs_hit10))
    std_diff_hit10 = float(np.std(test_diffs_hit10, ddof=1)) if len(test_diffs_hit10) > 1 else 0.0

    counts = Counter(best_triples_list)
    print("\n--- Selection Frequency of Best Triples across 20 Splits ---")
    for triple, cnt in counts.most_common():
        print(f"  Triple {triple}: {cnt}/20 ({cnt/20*100:.1f}%)")

    print(f"\nTest Difference (Best - Default) MRR: Mean = {mean_diff_mrr:+.5f}, Std = {std_diff_mrr:.5f}")
    print(f"Test Difference (Best - Default) Hit@10: Mean = {mean_diff_hit10:+.5f}, Std = {std_diff_hit10:.5f}")

    return {
        "objective": objective_name,
        "selection_freq": counts,
        "mean_diff_mrr": mean_diff_mrr,
        "std_diff_mrr": std_diff_mrr,
        "mean_diff_hit10": mean_diff_hit10,
        "std_diff_hit10": std_diff_hit10,
        "split_rows": split_rows,
    }

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
    triples = generate_weight_triples()

    res_p2b = evaluate_objective("p2b_dn", all_queries["P2b_dn_0.4"], matrices, triples)
    res_p1 = evaluate_objective("p1_dn", all_queries["P1_dn"], matrices, triples)

if __name__ == "__main__":
    main()
