"""
experiments/step18_tie_sensitivity.py

Requirement F2: Tie Sensitivity Analysis
- Evaluates Protocols P1, P4, and P4b for Methods M2 to M7.
- Evaluates Hit@10 and MRR under three tie rules:
    1. Average rank (primary): 1.0 + N_> + 0.5 * N_=
    2. Best rank (optimistic): 1.0 + N_>
    3. Worst rank (pessimistic): 1.0 + N_> + N_=
- Saves output to experiments/results/tie_sensitivity.csv.
- Inspects tie block sizes and rank distribution for M3 on P4 around ranks 10 to 20.
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

def compute_three_ranks(scores: np.ndarray, cand_indices: np.ndarray):
    """
    Computes average, best, and worst ranks for each candidate in cand_indices.
    Uses numpy sorting / searchsorted for exact O(M log M) performance.
    """
    cand_scores = scores[cand_indices]
    unique_vals, counts = np.unique(cand_scores, return_counts=True)
    # unique_vals is sorted in ascending order
    
    # Cumulative counts from highest to lowest
    # N_greater(u) = sum of counts of unique values > u
    n_unique = len(unique_vals)
    # suffix sum of counts
    suffix_sum = np.zeros(n_unique, dtype=np.int64)
    # suffix_sum[i] = sum(counts[i+1:])
    suffix_sum[:-1] = np.cumsum(counts[::-1])[:-1][::-1]
    
    # Map each unique value to (n_gt, n_eq)
    # n_eq = count - 1
    val_to_gt = {unique_vals[i]: suffix_sum[i] for i in range(n_unique)}
    val_to_eq = {unique_vals[i]: counts[i] - 1 for i in range(n_unique)}
    
    ranks_avg = {}
    ranks_best = {}
    ranks_worst = {}
    
    for idx in cand_indices:
        sc = scores[idx]
        n_gt = val_to_gt[sc]
        n_eq = val_to_eq[sc]
        ranks_best[idx] = 1.0 + float(n_gt)
        ranks_worst[idx] = 1.0 + float(n_gt + n_eq)
        ranks_avg[idx] = 1.0 + float(n_gt) + 0.5 * float(n_eq)
        
    return ranks_avg, ranks_best, ranks_worst

def run_tie_sensitivity():
    print("=" * 80)
    print("STEP F2: TIE SENSITIVITY ANALYSIS (PROTOCOLS P1, P4, P4b)")
    print("=" * 80)
    
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
    
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    pool_with_atc = sorted([cid for cid in pool_cids if len(cid_atc_l1[cid]) > 0])
    atc_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc], dtype=np.int32)
    
    methods = ["M2", "M3", "M4", "M5", "M6", "M7"]
    protocols = ["P1", "P4", "P4b"]
    
    results = []
    
    # Store per-query details for M3 on P4 to inspect rank distribution
    m3_p4_details = []
    
    for protocol in protocols:
        print(f"\nProcessing Protocol: {protocol}...")
        
        # Build query set
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
                    "q_cid": q_cid,
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
                    "q_cid": q_cid,
                    "q_idx": q_idx,
                    "cand_indices": cand_indices,
                    "rel_indices": rel_indices,
                })
                
        n_queries = len(queries)
        print(f"Protocol {protocol}: N = {n_queries} queries.")
        
        # Track metrics for each method and rule
        # rule -> method -> list of values
        hit10_data = {r: {m: [] for m in methods} for r in ["average", "best", "worst"]}
        hit20_data = {r: {m: [] for m in methods} for r in ["average", "best", "worst"]}
        mrr_data = {r: {m: [] for m in methods} for r in ["average", "best", "worst"]}
        
        for q in queries:
            q_idx = q["q_idx"]
            cands = q["cand_indices"]
            rels = q["rel_indices"]
            
            s = S[q_idx]
            t = T[q_idx]
            e = E[q_idx]
            m6 = 0.4 * s + 0.4 * t + 0.2 * e
            score_dict = {
                "M2": s,
                "M3": t,
                "M4": e,
                "M5": 0.5 * s + 0.5 * t,
                "M6": m6,
                "M7": m6 / penalties,
            }
            
            for m in methods:
                scores = score_dict[m]
                ranks_avg, ranks_best, ranks_worst = compute_three_ranks(scores, cands)
                
                # Minimum rank among relevant items
                r_avg_min = min(ranks_avg[c] for c in rels)
                r_best_min = min(ranks_best[c] for c in rels)
                r_worst_min = min(ranks_worst[c] for c in rels)
                
                hit10_data["average"][m].append(1.0 if r_avg_min <= 10.0 else 0.0)
                hit10_data["best"][m].append(1.0 if r_best_min <= 10.0 else 0.0)
                hit10_data["worst"][m].append(1.0 if r_worst_min <= 10.0 else 0.0)
                
                hit20_data["average"][m].append(1.0 if r_avg_min <= 20.0 else 0.0)
                hit20_data["best"][m].append(1.0 if r_best_min <= 20.0 else 0.0)
                hit20_data["worst"][m].append(1.0 if r_worst_min <= 20.0 else 0.0)
                
                mrr_data["average"][m].append(1.0 / r_avg_min)
                mrr_data["best"][m].append(1.0 / r_best_min)
                mrr_data["worst"][m].append(1.0 / r_worst_min)
                
                if protocol == "P4" and m == "M3":
                    # Record for analysis
                    m3_p4_details.append({
                        "q_cid": q["q_cid"],
                        "r_avg": r_avg_min,
                        "r_best": r_best_min,
                        "r_worst": r_worst_min,
                        "n_cands": len(cands),
                        "n_pos_t": sum(1 for c in cands if t[c] > 0.0),
                    })
                    
        for m in methods:
            for rule in ["average", "best", "worst"]:
                h10 = float(np.mean(hit10_data[rule][m]))
                h20 = float(np.mean(hit20_data[rule][m]))
                mrr = float(np.mean(mrr_data[rule][m]))
                results.append({
                    "protocol": protocol,
                    "method": m,
                    "tie_rule": rule,
                    "hit@10": round(h10, 4),
                    "hit@20": round(h20, 4),
                    "mrr": round(mrr, 4),
                    "n_queries": n_queries,
                })
                
    df_res = pd.DataFrame(results)
    out_csv = RESULTS_DIR / "tie_sensitivity.csv"
    df_res.to_csv(out_csv, index=False)
    print(f"\nSaved tie sensitivity results to: {out_csv}")
    
    # Analyze M3 on P4 rank jump
    df_m3 = pd.DataFrame(m3_p4_details)
    in_10_avg = (df_m3["r_avg"] <= 10.0).sum()
    in_20_avg = (df_m3["r_avg"] <= 20.0).sum()
    between_10_20 = ((df_m3["r_avg"] > 10.0) & (df_m3["r_avg"] <= 20.0)).sum()
    
    print("\n--- M3 on Protocol P4 Rank Distribution Audit ---")
    print(f"Total queries: {len(df_m3)}")
    print(f"Queries with r_avg <= 10.0: {in_10_avg} ({in_10_avg / len(df_m3):.4f})")
    print(f"Queries with r_avg <= 20.0: {in_20_avg} ({in_20_avg / len(df_m3):.4f})")
    print(f"Queries with 10.0 < r_avg <= 20.0: {between_10_20}")
    
    # Check what scores and tie block sizes those queries have
    sample_between = df_m3[(df_m3["r_avg"] > 10.0) & (df_m3["r_avg"] <= 20.0)]
    print("\nSample queries with 10 < r_avg <= 20:")
    print(sample_between.head(10))

if __name__ == "__main__":
    run_tie_sensitivity()
