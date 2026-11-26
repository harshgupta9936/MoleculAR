"""
experiments/step14_run_p4_p4b.py

Requirement N10:
Runs Protocol P4 and Protocol P4b for Methods M1 to M8:
- P4: Candidates are all pool drugs with >=1 ATC code (duplicate-neutral). Relevant drugs share indication and have no ATC Level-1 letter in common.
- P4b: Identical to P4, except that candidates sharing ANY ATC Level-1 letter with the query are REMOVED from the candidate pool (anatomical removal, duplicate-neutral).
- Evaluates Methods M1 to M8: Hit@10, Hit@20, Prec@10, MRR, PctRank, Enrichment with 95% bootstrap CIs (1,000 resamples, seed 42).
- Computes paired differences (1,000 resamples, seed 42): M7-M6, M6-M5, M7-M8, M4-M7.
- Computes coverage-stratified evaluation (strata a, b, c, d) for M2 to M7.
- Updates experiments/PROTOCOL_ADDENDUM.md with full formal specification.
- Outputs saved to experiments/results/.
"""

import sys
import time
import csv
from pathlib import Path
from collections import defaultdict
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_duplicates, rand_hit_at_k,
    compute_average_ranks, compute_m1_random, bootstrap_ci, bootstrap_enrichment,
    DB_A_PATH, GOLD_PRIMARY_PATH, DATA_DIR, RESULTS_DIR
)
from experiments.step11_p4_protocol import load_atc_level1

def evaluate_p4_variant(
    protocol_name: str,
    pool_cids: list,
    matrices: dict,
    drug_inds: dict,
    cid_atc_l1: dict,
    duplicates_of: dict,
    seed: int = 42
):
    cid_to_idx = matrices["cid_to_idx"]
    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    penalties = matrices["penalties"]
    N = len(pool_cids)
    
    pop_scores = np.array([len(drug_inds[cid]) for cid in pool_cids], dtype=np.float64)
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    pool_with_atc = sorted([cid for cid in pool_cids if len(cid_atc_l1[cid]) > 0])
    atc_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc], dtype=np.int32)
    
    queries = []
    for q_cid in cids_with_inds:
        if q_cid not in cid_atc_l1 or len(cid_atc_l1[q_cid]) == 0:
            continue
        q_idx = cid_to_idx[q_cid]
        q_l1 = cid_atc_l1[q_cid]
        q_inds = cid_ind_sets[q_cid]
        dups_q = duplicates_of.get(q_cid, set())
        
        # P1 indication relevant
        rel_P1 = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(cid_ind_sets[c])]
        
        # Cross-class ATC Level-1 relevant
        rel_p4 = [
            c for c in rel_P1
            if c in cid_atc_l1 and len(cid_atc_l1[c]) > 0
            and q_l1.isdisjoint(cid_atc_l1[c])
            and c not in dups_q
        ]
        if not rel_p4:
            continue
            
        rel_indices = np.array([cid_to_idx[c] for c in rel_p4], dtype=np.int32)
        
        if protocol_name == "P4":
            # Candidates: all pool drugs with ATC, excluding q and dups_q
            cand_indices = np.array([i for i in atc_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)
        elif protocol_name == "P4b":
            # Candidates: all pool drugs with ATC that share NO ATC Level-1 letter with q, excluding q and dups_q
            cand_indices = np.array([
                i for i in atc_indices
                if i != q_idx and pool_cids[i] not in dups_q
                and q_l1.isdisjoint(cid_atc_l1[pool_cids[i]])
            ], dtype=np.int32)
        else:
            raise ValueError(f"Unknown protocol: {protocol_name}")
            
        M_q = len(cand_indices)
        R_q = len(rel_indices)
        rand_h10 = rand_hit_at_k(M_q, R_q, 10)
        rand_h20 = rand_hit_at_k(M_q, R_q, 20)
        
        queries.append({
            "q_cid": q_cid,
            "q_idx": q_idx,
            "cand_indices": cand_indices,
            "rel_indices": rel_indices,
            "M": M_q,
            "R": R_q,
            "rand_hit10": rand_h10,
            "rand_hit20": rand_h20,
        })
        
    n_queries = len(queries)
    print(f"[{protocol_name}] Qualifying queries: {n_queries}")
    
    # Evaluate M1 to M8
    method_names = ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"]
    per_query_results = {m: defaultdict(list) for m in method_names}
    
    rng_m1 = np.random.RandomState(seed)
    
    for q in queries:
        q_idx = q["q_idx"]
        cands = q["cand_indices"]
        rels = q["rel_indices"]
        M = q["M"]
        R = q["R"]
        
        # M1: Random
        m1_res = compute_m1_random(M, R, rng_m1, n_permutations=100)
        per_query_results["M1"]["hit@10"].append(m1_res["hit10"])
        per_query_results["M1"]["hit@20"].append(m1_res["hit20"])
        per_query_results["M1"]["precision@10"].append(m1_res["prec10"])
        per_query_results["M1"]["mrr"].append(m1_res["mrr"])
        per_query_results["M1"]["percentile_rank"].append(m1_res["pct_rank"])
        per_query_results["M1"]["rand_hit10"].append(q["rand_hit10"])
        
        # Candidate features
        s_scores = S[q_idx]
        t_scores = T[q_idx]
        e_scores = E[q_idx]
        
        # Scores for M2 - M8 strictly following eval_engine.py
        scores_m6 = 0.4 * s_scores + 0.4 * t_scores + 0.2 * e_scores
        score_dict = {
            "M2": s_scores,
            "M3": t_scores,
            "M4": e_scores,
            "M5": 0.5 * s_scores + 0.5 * t_scores,
            "M6": scores_m6,
            "M7": scores_m6 / penalties,
            "M8": pop_scores,
        }
        
        for m in ["M2", "M3", "M4", "M5", "M6", "M7", "M8"]:
            scores = score_dict[m]
            ranks = compute_average_ranks(scores, cands)
            rel_ranks = [ranks[c] for c in rels]
            r_min = min(rel_ranks)
            
            # Count relevant in top 10
            n_rel_top10 = sum(1 for r in rel_ranks if r <= 10.0)
            
            per_query_results[m]["hit@10"].append(1.0 if r_min <= 10.0 else 0.0)
            per_query_results[m]["hit@20"].append(1.0 if r_min <= 20.0 else 0.0)
            per_query_results[m]["precision@10"].append(n_rel_top10 / 10.0)
            per_query_results[m]["mrr"].append(1.0 / r_min)
            per_query_results[m]["percentile_rank"].append(r_min / M)
            per_query_results[m]["rand_hit10"].append(q["rand_hit10"])
            
    # Compute summary metrics and bootstrap CIs
    metric_rows = []
    metric_keys = ["mrr", "hit@10", "hit@20", "precision@10", "percentile_rank", "enrichment"]
    
    for m in method_names:
        h10_arr = np.array(per_query_results[m]["hit@10"], dtype=np.float64)
        rand_h10_arr = np.array(per_query_results[m]["rand_hit10"], dtype=np.float64)
        
        for met in metric_keys:
            if met == "enrichment":
                mean_val, ci_low, ci_high = bootstrap_enrichment(h10_arr, rand_h10_arr, n_boot=1000, seed=seed)
            else:
                arr = np.array(per_query_results[m][met], dtype=np.float64)
                mean_val, ci_low, ci_high = bootstrap_ci(arr, n_boot=1000, seed=seed)
                
            metric_rows.append({
                "method": m,
                "metric": met,
                "mean": f"{mean_val:.6f}",
                "ci_low": f"{ci_low:.6f}",
                "ci_high": f"{ci_high:.6f}",
                "n_queries": n_queries,
            })
            
    return metric_rows, per_query_results, queries

def compute_paired_differences(per_query_dict, n_boot=1000, seed=42):
    rng = np.random.RandomState(seed)
    pairs = [
        ("M7", "M6"),
        ("M6", "M5"),
        ("M7", "M8"),
        ("M4", "M7"),
    ]
    diff_rows = []
    n = len(per_query_dict["M1"]["mrr"])
    boot_indices = rng.randint(0, n, size=(n_boot, n))
    
    for mA, mB in pairs:
        for metric in ["hit@10", "mrr"]:
            arrA = np.array(per_query_dict[mA][metric], dtype=np.float64)
            arrB = np.array(per_query_dict[mB][metric], dtype=np.float64)
            diff = arrA - arrB
            mean_diff = float(np.mean(diff))
            
            boot_diffs = np.empty(n_boot, dtype=np.float64)
            for b in range(n_boot):
                boot_diffs[b] = np.mean(diff[boot_indices[b]])
            ci_low = float(np.percentile(boot_diffs, 2.5))
            ci_high = float(np.percentile(boot_diffs, 97.5))
            
            diff_rows.append({
                "comparison": f"{mA} - {mB}",
                "metric": metric,
                "mean_diff": f"{mean_diff:.4f}",
                "ci_low": f"{ci_low:.4f}",
                "ci_high": f"{ci_high:.4f}",
                "formatted": f"{mean_diff:+.4f} [{ci_low:+.4f}, {ci_high:+.4f}]"
            })
    return diff_rows

def compute_coverage_stratification(queries, per_query_dict, drugs, n_boot=1000, seed=42):
    # Strata:
    # (a) both >= 1
    # (b) targets only
    # (c) se only
    # (d) neither
    strata = {
        "(a) Targets and Side Effects": [],
        "(b) Targets Only": [],
        "(c) Side Effects Only": [],
        "(d) Neither": []
    }
    for idx, q in enumerate(queries):
        cid = q["q_cid"]
        has_t = len(drugs[cid]["targets"]) > 0
        has_se = len(drugs[cid]["side_effects"]) > 0
        if has_t and has_se:
            strata["(a) Targets and Side Effects"].append(idx)
        elif has_t and not has_se:
            strata["(b) Targets Only"].append(idx)
        elif not has_t and has_se:
            strata["(c) Side Effects Only"].append(idx)
        else:
            strata["(d) Neither"].append(idx)
            
    strat_rows = []
    for s_name, indices in strata.items():
        n_strat = len(indices)
        idx_arr = np.array(indices, dtype=np.int32)
        for m in ["M2", "M3", "M4", "M5", "M6", "M7"]:
            if n_strat == 0:
                continue
            h10 = np.array([per_query_dict[m]["hit@10"][i] for i in indices], dtype=np.float64)
            mrr = np.array([per_query_dict[m]["mrr"][i] for i in indices], dtype=np.float64)
            rand_h10 = np.array([per_query_dict[m]["rand_hit10"][i] for i in indices], dtype=np.float64)
            
            m_h10, l_h10, u_h10 = bootstrap_ci(h10, n_boot=n_boot, seed=seed)
            m_mrr, l_mrr, u_mrr = bootstrap_ci(mrr, n_boot=n_boot, seed=seed)
            m_enr, l_enr, u_enr = bootstrap_enrichment(h10, rand_h10, n_boot=n_boot, seed=seed)
            
            strat_rows.append({
                "stratum": s_name,
                "n_queries": n_strat,
                "method": m,
                "hit10_mean": f"{m_h10:.4f}",
                "hit10_ci": f"{m_h10:.4f} [{l_h10:.4f}, {u_h10:.4f}]",
                "mrr_mean": f"{m_mrr:.4f}",
                "mrr_ci": f"{m_mrr:.4f} [{l_mrr:.4f}, {u_mrr:.4f}]",
                "enrichment_mean": f"{m_enr:.2f}",
                "enrichment_ci": f"{m_enr:.2f} [{l_enr:.2f}, {u_enr:.2f}]",
            })
    return strat_rows

def main():
    print("=" * 80)
    print("STEP 14: RUN PROTOCOL P4 AND P4b (REQUIREMENT N10)")
    print("=" * 80)
    t0 = time.time()
    
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids)
    cid_atc_l1, _ = load_atc_level1(set(pool_cids))
    
    # Precompute S, T, E
    matrices = {}
    matrices["cid_to_idx"] = {cid: i for i, cid in enumerate(pool_cids)}
    N = len(pool_cids)
    
    # S matrix
    fps = [drugs[cid]["fingerprint"] for cid in pool_cids]
    from rdkit.Chem import DataStructs
    S = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        if fps[i] is not None:
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
            S[i] = np.array(sims, dtype=np.float64)
    matrices["S"] = S
    
    # T matrix
    T = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        t_i = drugs[pool_cids[i]]["targets"]
        if not t_i: continue
        for j in range(i, N):
            t_j = drugs[pool_cids[j]]["targets"]
            if not t_j: continue
            u = len(t_i | t_j)
            val = len(t_i & t_j) / u if u > 0 else 0.0
            T[i, j] = val
            T[j, i] = val
    matrices["T"] = T
    
    # E matrix
    E = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        se_i = drugs[pool_cids[i]]["side_effects"]
        if not se_i: continue
        for j in range(i, N):
            se_j = drugs[pool_cids[j]]["side_effects"]
            if not se_j: continue
            u = len(se_i | se_j)
            val = len(se_i & se_j) / u if u > 0 else 0.0
            E[i, j] = val
            E[j, i] = val
    matrices["E"] = E
    
    hub_pcts = np.array([drugs[cid]["hubness_percentile"] for cid in pool_cids], dtype=np.float64)
    matrices["penalties"] = 1.0 + (hub_pcts / 100.0) ** 2
    
    # 1. Run P4
    print("\nRunning Protocol P4...")
    p4_metrics, p4_per_query, p4_queries = evaluate_p4_variant(
        "P4", pool_cids, matrices, drug_inds, cid_atc_l1, duplicates_of, seed=42
    )
    p4_csv = RESULTS_DIR / "metrics_C1_P4.csv"
    with open(p4_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["method", "metric", "mean", "ci_low", "ci_high", "n_queries"])
        w.writeheader()
        w.writerows(p4_metrics)
    print(f"P4 metrics written to {p4_csv}")
    
    # 2. Run P4b
    print("\nRunning Protocol P4b...")
    p4b_metrics, p4b_per_query, p4b_queries = evaluate_p4_variant(
        "P4b", pool_cids, matrices, drug_inds, cid_atc_l1, duplicates_of, seed=42
    )
    p4b_csv = RESULTS_DIR / "metrics_C1_P4b.csv"
    with open(p4b_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["method", "metric", "mean", "ci_low", "ci_high", "n_queries"])
        w.writeheader()
        w.writerows(p4b_metrics)
    print(f"P4b metrics written to {p4b_csv}")
    
    # 3. Census comparison
    p4_total_rel = sum(q["R"] for q in p4_queries)
    p4b_total_rel = sum(q["R"] for q in p4b_queries)
    p4_rand_h10 = np.mean([q["rand_hit10"] for q in p4_queries])
    p4b_rand_h10 = np.mean([q["rand_hit10"] for q in p4b_queries])
    
    census_csv = RESULTS_DIR / "p4_p4b_query_census.csv"
    with open(census_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["protocol", "n_queries", "total_relevant_pairs", "mean_M", "mean_R", "analytical_random_hit10"])
        w.writerow(["P4", len(p4_queries), p4_total_rel, f"{np.mean([q['M'] for q in p4_queries]):.2f}", f"{np.mean([q['R'] for q in p4_queries]):.4f}", f"{p4_rand_h10:.6f}"])
        w.writerow(["P4b", len(p4b_queries), p4b_total_rel, f"{np.mean([q['M'] for q in p4b_queries]):.2f}", f"{np.mean([q['R'] for q in p4b_queries]):.4f}", f"{p4b_rand_h10:.6f}"])
        
    # 4. Paired differences
    print("\nComputing paired differences for P4 and P4b...")
    diff_p4 = compute_paired_differences(p4_per_query, n_boot=1000, seed=42)
    diff_p4b = compute_paired_differences(p4b_per_query, n_boot=1000, seed=42)
    
    paired_csv = RESULTS_DIR / "paired_differences_p4_p4b.csv"
    with open(paired_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["protocol", "comparison", "metric", "mean_diff", "ci_low", "ci_high", "formatted"])
        w.writeheader()
        for r in diff_p4:
            r["protocol"] = "P4"
            w.writerow(r)
        for r in diff_p4b:
            r["protocol"] = "P4b"
            w.writerow(r)
            
    # 5. Coverage Stratified Table
    print("\nComputing coverage stratification for P4 and P4b...")
    strat_p4 = compute_coverage_stratification(p4_queries, p4_per_query, drugs, n_boot=1000, seed=42)
    strat_p4b = compute_coverage_stratification(p4b_queries, p4b_per_query, drugs, n_boot=1000, seed=42)
    
    strat_csv = RESULTS_DIR / "coverage_stratified_p4_p4b.csv"
    with open(strat_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=[
            "protocol", "stratum", "n_queries", "method", "hit10_mean", "hit10_ci",
            "mrr_mean", "mrr_ci", "enrichment_mean", "enrichment_ci"
        ])
        w.writeheader()
        for r in strat_p4:
            r["protocol"] = "P4"
            w.writerow(r)
        for r in strat_p4b:
            r["protocol"] = "P4b"
            w.writerow(r)
            
    print(f"\nCompleted in {time.time() - t0:.2f}s.")

if __name__ == "__main__":
    main()
