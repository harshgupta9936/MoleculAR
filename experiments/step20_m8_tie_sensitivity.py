"""
experiments/step20_m8_tie_sensitivity.py

Requirement F6: M8 Tie Sensitivity & Paired Differences M7 - M8
- Extends tie sensitivity analysis to Method M8 (Popularity) for Protocols P1, P4, and P4b:
    * Best rank (optimistic): 1.0 + N_>
    * Average rank (primary): 1.0 + N_> + 0.5 * N_=
    * Worst rank (pessimistic): 1.0 + N_> + N_=
    * Metrics: Hit@10 and MRR
- Computes paired differences M7 minus M8 on Protocols P4 and P4b under each of the three tie rules:
    * Metric difference: delta_i = Metric(M7)_i - Metric(M8)_i
    * 1,000 bootstrap resamples with seed 42 for 95% confidence intervals
- Saves output to experiments/results/tie_sensitivity_m8.csv
- Label: SENSITIVITY ONLY (Primary tables unchanged)
"""

import sys
import csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_duplicates,
    precompute_matrices, DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)
from experiments.step11_p4_protocol import load_atc_level1
from experiments.step18_tie_sensitivity import compute_three_ranks

def bootstrap_paired_diff(diff_arr: np.ndarray, n_boot: int = 1000, seed: int = 42):
    rng = np.random.RandomState(seed)
    n = len(diff_arr)
    boot_means = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        idx = rng.randint(0, n, size=n)
        boot_means[i] = np.mean(diff_arr[idx])
    mean_val = float(np.mean(diff_arr))
    ci_low = float(np.percentile(boot_means, 2.5))
    ci_high = float(np.percentile(boot_means, 97.5))
    formatted = f"{mean_val:+.4f} [{ci_low:+.4f}, {ci_high:+.4f}]"
    return mean_val, ci_low, ci_high, formatted

def run_m8_tie_sensitivity():
    print("=" * 85)
    print("STEP F6: M8 TIE SENSITIVITY & PAIRED DIFFERENCES (M7 - M8)")
    print("=" * 85)

    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    matrices = precompute_matrices(drugs, pool_cids)
    cid_to_idx = matrices["cid_to_idx"]
    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    penalties = matrices["penalties"]
    all_pool_indices = np.arange(len(pool_cids), dtype=np.int32)

    cid_atc_l1, _ = load_atc_level1(set(pool_cids))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids)

    pop_scores = np.array([len(drug_inds[c]) for c in pool_cids], dtype=np.float64)
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    pool_with_atc = sorted([cid for cid in pool_cids if len(cid_atc_l1[cid]) > 0])
    atc_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc], dtype=np.int32)

    protocols = ["P1", "P4", "P4b"]
    rules = ["average", "best", "worst"]

    # Store per-query scores: protocol -> rule -> method -> list of per-query metric
    per_query_metrics = {
        p: {r: {"M7": {"hit@10": [], "mrr": []}, "M8": {"hit@10": [], "mrr": []}} for r in rules}
        for p in protocols
    }

    n_queries_map = {}

    for protocol in protocols:
        print(f"Processing {protocol}...")
        queries = []
        if protocol == "P1":
            for q_cid in cids_with_inds:
                q_idx = cid_to_idx[q_cid]
                q_inds = cid_ind_sets[q_cid]
                rels = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(cid_ind_sets[c])]
                if not rels:
                    continue
                cand_indices = np.array([i for i in all_pool_indices if i != q_idx], dtype=np.int32)
                rel_indices = np.array([cid_to_idx[c] for c in rels], dtype=np.int32)
                queries.append({
                    "q_idx": q_idx,
                    "cand_indices": cand_indices,
                    "rel_indices": rel_indices,
                })
        else: # P4 and P4b
            for q_cid in cids_with_inds:
                if q_cid not in cid_atc_l1 or len(cid_atc_l1[q_cid]) == 0:
                    continue
                q_idx = cid_to_idx[q_cid]
                q_l1 = cid_atc_l1[q_cid]
                q_inds = cid_ind_sets[q_cid]
                dups_q = duplicates_of.get(q_cid, set())
                rels = [
                    c for c in cids_with_inds
                    if c != q_cid and c in cid_atc_l1 and len(cid_atc_l1[c]) > 0
                    and q_l1.isdisjoint(cid_atc_l1[c]) and c not in dups_q
                    and not q_inds.isdisjoint(cid_ind_sets[c])
                ]
                if not rels:
                    continue
                rel_indices = np.array([cid_to_idx[c] for c in rels], dtype=np.int32)
                if protocol == "P4":
                    cand_indices = np.array([i for i in atc_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)
                else: # P4b
                    cand_indices = np.array([
                        i for i in atc_indices
                        if i != q_idx and pool_cids[i] not in dups_q
                        and q_l1.isdisjoint(cid_atc_l1[pool_cids[i]])
                    ], dtype=np.int32)
                queries.append({
                    "q_idx": q_idx,
                    "cand_indices": cand_indices,
                    "rel_indices": rel_indices,
                })

        n_queries_map[protocol] = len(queries)

        for q in queries:
            q_idx = q["q_idx"]
            cands = q["cand_indices"]
            rels = q["rel_indices"]

            s = S[q_idx]
            t = T[q_idx]
            e = E[q_idx]
            m6 = 0.4 * s + 0.4 * t + 0.2 * e
            m7 = m6 / penalties
            m8 = pop_scores

            # Ranks for M7
            r_avg_m7, r_best_m7, r_worst_m7 = compute_three_ranks(m7, cands)
            # Ranks for M8
            r_avg_m8, r_best_m8, r_worst_m8 = compute_three_ranks(m8, cands)

            for rule, r_m7_map, r_m8_map in [
                ("average", r_avg_m7, r_avg_m8),
                ("best", r_best_m7, r_best_m8),
                ("worst", r_worst_m7, r_worst_m8),
            ]:
                min_r_m7 = min(r_m7_map[c] for c in rels)
                min_r_m8 = min(r_m8_map[c] for c in rels)

                per_query_metrics[protocol][rule]["M7"]["hit@10"].append(1.0 if min_r_m7 <= 10.0 else 0.0)
                per_query_metrics[protocol][rule]["M7"]["mrr"].append(1.0 / min_r_m7)

                per_query_metrics[protocol][rule]["M8"]["hit@10"].append(1.0 if min_r_m8 <= 10.0 else 0.0)
                per_query_metrics[protocol][rule]["M8"]["mrr"].append(1.0 / min_r_m8)

    # Compile results table
    rows = []

    # 1. M8 Absolute Performance under the 3 tie rules
    for protocol in protocols:
        nq = n_queries_map[protocol]
        for rule in rules:
            h10 = float(np.mean(per_query_metrics[protocol][rule]["M8"]["hit@10"]))
            mrr = float(np.mean(per_query_metrics[protocol][rule]["M8"]["mrr"]))
            rows.append({
                "category": "m8_performance",
                "protocol": protocol,
                "comparison_or_method": "M8",
                "tie_rule": rule,
                "metric": "hit@10",
                "value_or_diff": round(h10, 4),
                "ci_low": np.nan,
                "ci_high": np.nan,
                "formatted": f"{h10:.4f}",
                "n_queries": nq,
            })
            rows.append({
                "category": "m8_performance",
                "protocol": protocol,
                "comparison_or_method": "M8",
                "tie_rule": rule,
                "metric": "mrr",
                "value_or_diff": round(mrr, 4),
                "ci_low": np.nan,
                "ci_high": np.nan,
                "formatted": f"{mrr:.4f}",
                "n_queries": nq,
            })

    # 2. Paired differences M7 minus M8 on P4 and P4b under each of the 3 tie rules
    for protocol in ["P4", "P4b"]:
        nq = n_queries_map[protocol]
        for rule in rules:
            for met in ["hit@10", "mrr"]:
                arr_m7 = np.array(per_query_metrics[protocol][rule]["M7"][met], dtype=np.float64)
                arr_m8 = np.array(per_query_metrics[protocol][rule]["M8"][met], dtype=np.float64)
                diff_arr = arr_m7 - arr_m8
                m_diff, ci_l, ci_h, fmt = bootstrap_paired_diff(diff_arr, n_boot=1000, seed=42)
                rows.append({
                    "category": "paired_difference_m7_minus_m8",
                    "protocol": protocol,
                    "comparison_or_method": "M7 - M8",
                    "tie_rule": rule,
                    "metric": met,
                    "value_or_diff": round(m_diff, 4),
                    "ci_low": round(ci_l, 4),
                    "ci_high": round(ci_h, 4),
                    "formatted": fmt,
                    "n_queries": nq,
                })

    df_out = pd.DataFrame(rows)
    out_csv = RESULTS_DIR / "tie_sensitivity_m8.csv"
    df_out.to_csv(out_csv, index=False)
    print(f"\nSaved M8 tie sensitivity to: {out_csv}")
    print("\n--- Output Summary ---")
    print(df_out.to_string())

if __name__ == "__main__":
    run_m8_tie_sensitivity()
