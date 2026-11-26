"""
experiments/step15_coverage_fair_nulls.py

Requirement N11:
Coverage-Fair Nulls (Shuffle Control within Coverage Strata):
- Re-runs the shuffle control with the same permutations as N7 (20 runs, seed 42).
- M3 evaluated within coverage strata (a) and (b), plus combined targets-present.
- M4 evaluated within coverage strata (a) and (c), plus combined side-effects-present.
- Protocols: P1, P2a(0.4), P3-dn, and P4.
- Reports Real vs Null Mean and Range [Min, Max] for Hit@10 and MRR.
- Does NOT report enrichment for single-signal methods in this table.
- Output saved to experiments/results/coverage_fair_nulls.csv.
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
    load_database, load_gold_pairs, load_duplicates, load_atc_l3,
    compute_average_ranks, DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)
from experiments.step10_shuffle_control import compute_jaccard_matrix
from experiments.step11_p4_protocol import load_atc_level1

def evaluate_method_per_query(sim_matrix: np.ndarray, queries: list) -> tuple:
    hits = np.empty(len(queries), dtype=np.float64)
    mrrs = np.empty(len(queries), dtype=np.float64)
    for idx, (q_idx, cand_indices, rel_indices) in enumerate(queries):
        ranks_dict = compute_average_ranks(sim_matrix[q_idx], cand_indices)
        rel_ranks = [ranks_dict[c] for c in rel_indices]
        r_min = min(rel_ranks)
        hits[idx] = 1.0 if r_min <= 10.0 else 0.0
        mrrs[idx] = 1.0 / r_min
    return hits, mrrs

def main():
    print("=" * 80)
    print("STEP 15: COVERAGE-FAIR NULLS (REQUIREMENT N11)")
    print("=" * 80)
    t0 = time.time()
    
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    N = len(pool_cids)
    cid_to_idx = {cid: i for i, cid in enumerate(pool_cids)}
    
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    cid_atc_l3 = load_atc_l3(set(pool_cids))
    cid_atc_l1, _ = load_atc_level1(set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids)
    
    # Precompute unpermuted T and E
    T_real = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        t_i = drugs[pool_cids[i]]["targets"]
        if not t_i: continue
        for j in range(i, N):
            t_j = drugs[pool_cids[j]]["targets"]
            if not t_j: continue
            u = len(t_i | t_j)
            val = len(t_i & t_j) / u if u > 0 else 0.0
            T_real[i, j] = val
            T_real[j, i] = val
            
    E_real = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        se_i = drugs[pool_cids[i]]["side_effects"]
        if not se_i: continue
        for j in range(i, N):
            se_j = drugs[pool_cids[j]]["side_effects"]
            if not se_j: continue
            u = len(se_i | se_j)
            val = len(se_i & se_j) / u if u > 0 else 0.0
            E_real[i, j] = val
            E_real[j, i] = val

    # S matrix for P2a
    fps = [drugs[cid]["fingerprint"] for cid in pool_cids]
    from rdkit.Chem import DataStructs
    S = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        if fps[i] is not None:
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
            S[i] = np.array(sims, dtype=np.float64)

    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    pool_with_atc_l3 = {cid for cid in pool_cids if len(cid_atc_l3[cid]) > 0}
    atc_l3_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc_l3], dtype=np.int32)
    pool_with_atc_l1 = {cid for cid in pool_cids if len(cid_atc_l1[cid]) > 0}
    atc_l1_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc_l1], dtype=np.int32)

    protocols = ["P1", "P2a_0.4", "P3_dn", "P4"]
    # We will build query lists for each protocol
    # and classify each query into its coverage stratum
    proto_queries = {p: [] for p in protocols}

    for q_cid in cids_with_inds:
        q_idx = cid_to_idx[q_cid]
        q_inds = cid_ind_sets[q_cid]
        dups_q = duplicates_of.get(q_cid, set())
        has_t = len(drugs[q_cid]["targets"]) > 0
        has_se = len(drugs[q_cid]["side_effects"]) > 0
        
        # Determine stratum
        if has_t and has_se:
            strat = "(a) Targets and Side Effects"
        elif has_t and not has_se:
            strat = "(b) Targets Only"
        elif not has_t and has_se:
            strat = "(c) Side Effects Only"
        else:
            strat = "(d) Neither"
            
        rel_P1 = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(cid_ind_sets[c])]
        if not rel_P1:
            continue
            
        # P1
        cand_p1 = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)
        rel_p1_indices = np.array([cid_to_idx[c] for c in rel_P1], dtype=np.int32)
        proto_queries["P1"].append((q_idx, cand_p1, rel_p1_indices, strat))
        
        # P2a_0.4
        rel_p2a = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4]
        if rel_p2a:
            cand_p2a = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)
            rel_p2a_indices = np.array([cid_to_idx[c] for c in rel_p2a], dtype=np.int32)
            proto_queries["P2a_0.4"].append((q_idx, cand_p2a, rel_p2a_indices, strat))
            
        # P3_dn
        if q_cid in pool_with_atc_l3:
            q_atc = cid_atc_l3[q_cid]
            rel_p3 = [c for c in rel_P1 if c in pool_with_atc_l3 and q_atc.isdisjoint(cid_atc_l3[c]) and c not in dups_q]
            if rel_p3:
                cand_p3 = np.array([i for i in atc_l3_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)
                rel_p3_indices = np.array([cid_to_idx[c] for c in rel_p3], dtype=np.int32)
                proto_queries["P3_dn"].append((q_idx, cand_p3, rel_p3_indices, strat))
                
        # P4
        if q_cid in pool_with_atc_l1:
            q_l1 = cid_atc_l1[q_cid]
            rel_p4 = [c for c in rel_P1 if c in pool_with_atc_l1 and q_l1.isdisjoint(cid_atc_l1[c]) and c not in dups_q]
            if rel_p4:
                cand_p4 = np.array([i for i in atc_l1_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)
                rel_p4_indices = np.array([cid_to_idx[c] for c in rel_p4], dtype=np.int32)
                proto_queries["P4"].append((q_idx, cand_p4, rel_p4_indices, strat))

    for p in protocols:
        print(f"Loaded {p}: {len(proto_queries[p])} queries")

    # Group queries by protocol and stratum
    proto_strata_queries = {}
    for p in protocols:
        proto_strata_queries[p] = defaultdict(list)
        for item in proto_queries[p]:
            strat = item[3]
            proto_strata_queries[p][strat].append((item[0], item[1], item[2]))

    # Target drugs & SE drugs for shuffling (identical to N7)
    target_drugs_idx = [i for i, cid in enumerate(pool_cids) if len(drugs[cid]["targets"]) > 0]
    n_target_drugs = len(target_drugs_idx)
    orig_target_sets = [drugs[cid]["targets"] for cid in pool_cids]

    se_drugs_idx = [i for i, cid in enumerate(pool_cids) if len(drugs[cid]["side_effects"]) > 0]
    n_se_drugs = len(se_drugs_idx)
    orig_se_sets = [drugs[cid]["side_effects"] for cid in pool_cids]

    N_PERMS = 20
    rng_t = np.random.RandomState(42)
    rng_se = np.random.RandomState(42)

    # We evaluate:
    # For M3: stratum (a), stratum (b), and combined (a)+(b)
    # For M4: stratum (a), stratum (c), and combined (a)+(c)
    m3_eval_groups = ["(a) Targets and Side Effects", "(b) Targets Only", "(a)+(b) Combined (Targets Present)"]
    m4_eval_groups = ["(a) Targets and Side Effects", "(c) Side Effects Only", "(a)+(c) Combined (Side Effects Present)"]

    # 1. Real evaluations
    real_results = {}
    for p in protocols:
        # M3 real
        for grp in m3_eval_groups:
            if grp == "(a)+(b) Combined (Targets Present)":
                q_list = proto_strata_queries[p]["(a) Targets and Side Effects"] + proto_strata_queries[p]["(b) Targets Only"]
            else:
                q_list = proto_strata_queries[p][grp]
            if q_list:
                hits, mrrs = evaluate_method_per_query(T_real, q_list)
                real_results[(p, "M3", grp)] = (float(np.mean(hits)), float(np.mean(mrrs)), len(q_list))
            else:
                real_results[(p, "M3", grp)] = (0.0, 0.0, 0)
                
        # M4 real
        for grp in m4_eval_groups:
            if grp == "(a)+(c) Combined (Side Effects Present)":
                q_list = proto_strata_queries[p]["(a) Targets and Side Effects"] + proto_strata_queries[p]["(c) Side Effects Only"]
            else:
                q_list = proto_strata_queries[p][grp]
            if q_list:
                hits, mrrs = evaluate_method_per_query(E_real, q_list)
                real_results[(p, "M4", grp)] = (float(np.mean(hits)), float(np.mean(mrrs)), len(q_list))
            else:
                real_results[(p, "M4", grp)] = (0.0, 0.0, 0)

    # 2. Permutation runs
    print("\nRunning M3 Target Coverage-Fair Permutations (20 runs, seed 42)...")
    m3_null_runs = {k: {"hit10": [], "mrr": []} for k in real_results if k[1] == "M3"}
    for perm_idx in range(N_PERMS):
        perm_sub = rng_t.permutation(n_target_drugs)
        shuffled_targets = [set() for _ in range(N)]
        for i_sub, j_sub in enumerate(perm_sub):
            shuffled_targets[target_drugs_idx[i_sub]] = orig_target_sets[target_drugs_idx[j_sub]]
        T_perm = compute_jaccard_matrix(shuffled_targets, N)
        for p in protocols:
            for grp in m3_eval_groups:
                if grp == "(a)+(b) Combined (Targets Present)":
                    q_list = proto_strata_queries[p]["(a) Targets and Side Effects"] + proto_strata_queries[p]["(b) Targets Only"]
                else:
                    q_list = proto_strata_queries[p][grp]
                if q_list:
                    hits, mrrs = evaluate_method_per_query(T_perm, q_list)
                    m3_null_runs[(p, "M3", grp)]["hit10"].append(float(np.mean(hits)))
                    m3_null_runs[(p, "M3", grp)]["mrr"].append(float(np.mean(mrrs)))
        print(f"  M3 Permutation {perm_idx+1}/{N_PERMS} complete.")

    print("\nRunning M4 Side-Effect Coverage-Fair Permutations (20 runs, seed 42)...")
    m4_null_runs = {k: {"hit10": [], "mrr": []} for k in real_results if k[1] == "M4"}
    for perm_idx in range(N_PERMS):
        perm_sub = rng_se.permutation(n_se_drugs)
        shuffled_se = [set() for _ in range(N)]
        for i_sub, j_sub in enumerate(perm_sub):
            shuffled_se[se_drugs_idx[i_sub]] = orig_se_sets[se_drugs_idx[j_sub]]
        E_perm = compute_jaccard_matrix(shuffled_se, N)
        for p in protocols:
            for grp in m4_eval_groups:
                if grp == "(a)+(c) Combined (Side Effects Present)":
                    q_list = proto_strata_queries[p]["(a) Targets and Side Effects"] + proto_strata_queries[p]["(c) Side Effects Only"]
                else:
                    q_list = proto_strata_queries[p][grp]
                if q_list:
                    hits, mrrs = evaluate_method_per_query(E_perm, q_list)
                    m4_null_runs[(p, "M4", grp)]["hit10"].append(float(np.mean(hits)))
                    m4_null_runs[(p, "M4", grp)]["mrr"].append(float(np.mean(mrrs)))
        print(f"  M4 Permutation {perm_idx+1}/{N_PERMS} complete.")

    # 3. Save to experiments/results/coverage_fair_nulls.csv
    out_csv = RESULTS_DIR / "coverage_fair_nulls.csv"
    rows = []
    
    all_runs = {}
    all_runs.update(m3_null_runs)
    all_runs.update(m4_null_runs)
    
    for (p, m, grp), null_data in sorted(all_runs.items()):
        real_h10, real_mrr, n_q = real_results[(p, m, grp)]
        if n_q == 0:
            continue
        null_h10 = null_data["hit10"]
        null_mrr = null_data["mrr"]
        
        m_null_h10 = float(np.mean(null_h10))
        min_null_h10 = float(np.min(null_h10))
        max_null_h10 = float(np.max(null_h10))
        
        m_null_mrr = float(np.mean(null_mrr))
        min_null_mrr = float(np.min(null_mrr))
        max_null_mrr = float(np.max(null_mrr))
        
        rows.append({
            "protocol": p,
            "method": m,
            "stratum": grp,
            "n_queries": n_q,
            "real_hit10": f"{real_h10:.4f}",
            "null_mean_hit10": f"{m_null_h10:.4f}",
            "null_range_hit10": f"[{min_null_h10:.4f}, {max_null_h10:.4f}]",
            "real_mrr": f"{real_mrr:.4f}",
            "null_mean_mrr": f"{m_null_mrr:.4f}",
            "null_range_mrr": f"[{min_null_mrr:.4f}, {max_null_mrr:.4f}]",
            "n_permutations": N_PERMS,
            "random_seed": 42
        })
        
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "protocol", "method", "stratum", "n_queries",
            "real_hit10", "null_mean_hit10", "null_range_hit10",
            "real_mrr", "null_mean_mrr", "null_range_mrr",
            "n_permutations", "random_seed"
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
        
    print(f"\nCoverage-fair nulls complete in {time.time() - t0:.2f}s.")
    print(f"Results saved to {out_csv}")

if __name__ == "__main__":
    main()
