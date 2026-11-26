"""
experiments/step10_shuffle_control.py

Requirement N7:
Shuffle Control (sanity check) for P1, P2a(0.4), and P3-dn, M3 and M4 only:
- Randomly permute the target sets across drugs that have targets (20 permutations, seed 42)
- Separately permute the side-effect sets across drugs that have side effects (20 permutations, seed 42)
- Recompute hit@10, MRR, and enrichment with the same tie rule (compute_average_ranks)
- Compute mean and range [min, max] over the 20 runs
- Compare next to real values from on-disk CSVs
- Output saved to experiments/results/shuffle_control_results.csv
"""

import sys
import time
import csv
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_duplicates, load_atc_l3,
    precompute_matrices, compute_average_ranks, rand_hit_at_k,
    DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

def compute_jaccard_matrix(feature_sets: list, N: int) -> np.ndarray:
    mat = np.zeros((N, N), dtype=np.float64)
    # only compute for non-empty
    non_empty_indices = [i for i in range(N) if len(feature_sets[i]) > 0]
    for i_idx, i in enumerate(non_empty_indices):
        set_i = feature_sets[i]
        for j in non_empty_indices[i_idx:]:
            set_j = feature_sets[j]
            u = len(set_i | set_j)
            if u > 0:
                val = len(set_i & set_j) / u
                mat[i, j] = val
                mat[j, i] = val
    return mat

def evaluate_method_on_queries(sim_matrix: np.ndarray, queries: list) -> dict:
    hits = []
    mrrs = []
    rand_hits = []
    
    for q_idx, cand_indices, rel_indices, rand_h10 in queries:
        cand_scores = sim_matrix[q_idx, cand_indices]
        ranks_dict = compute_average_ranks(sim_matrix[q_idx], cand_indices)
        rel_ranks = [ranks_dict[c] for c in rel_indices]
        r_min = min(rel_ranks)
        
        hits.append(1.0 if r_min <= 10.0 else 0.0)
        mrrs.append(1.0 / r_min)
        rand_hits.append(rand_h10)
        
    m_hit = float(np.mean(hits))
    m_mrr = float(np.mean(mrrs))
    m_rand = float(np.mean(rand_hits))
    enrichment = m_hit / m_rand if m_rand > 0 else 1.0
    
    return {
        "hit10": m_hit,
        "mrr": m_mrr,
        "enrichment": enrichment,
    }

def main():
    print("=" * 80)
    print("STEP 10: SHUFFLE CONTROL FOR M3 AND M4 (REQUIREMENT N7)")
    print("=" * 80)
    t0 = time.time()
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Load Data
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    N = len(pool_cids)
    cid_to_idx = {cid: i for i, cid in enumerate(pool_cids)}
    
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    cid_atc_l3 = load_atc_l3(set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids)
    
    # Precompute standard S matrix for P2a
    fps = [drugs[cid]["fingerprint"] for cid in pool_cids]
    from rdkit.Chem import DataStructs
    S = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        if fps[i] is not None:
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
            S[i] = np.array(sims, dtype=np.float64)
            
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    pool_with_atc = {cid for cid in pool_cids if len(cid_atc_l3[cid]) > 0}
    atc_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc], dtype=np.int32)
    
    # Build query definitions for P1, P2a_0.4, P3_dn
    protocols = ["P1", "P2a_0.4", "P3_dn"]
    query_specs = {p: [] for p in protocols}
    
    for q_cid in cids_with_inds:
        q_idx = cid_to_idx[q_cid]
        q_inds = cid_ind_sets[q_cid]
        rel_P1 = [c_cid for c_cid in cids_with_inds if c_cid != q_cid and not q_inds.isdisjoint(cid_ind_sets[c_cid])]
        if not rel_P1:
            continue
            
        # P1
        cand_p1 = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)
        rel_p1_indices = [cid_to_idx[c] for c in rel_P1]
        rand_p1 = rand_hit_at_k(len(cand_p1), len(rel_p1_indices), 10)
        query_specs["P1"].append((q_idx, cand_p1, rel_p1_indices, rand_p1))
        
        # P2a_0.4
        rel_p2a = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4]
        if rel_p2a:
            cand_p2a = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)
            rel_p2a_indices = [cid_to_idx[c] for c in rel_p2a]
            rand_p2a = rand_hit_at_k(len(cand_p2a), len(rel_p2a_indices), 10)
            query_specs["P2a_0.4"].append((q_idx, cand_p2a, rel_p2a_indices, rand_p2a))
            
        # P3_dn
        if q_cid in pool_with_atc:
            q_atc = cid_atc_l3[q_cid]
            dups_q = duplicates_of.get(q_cid, set())
            rel_p3 = [c for c in rel_P1 if c in pool_with_atc and q_atc.isdisjoint(cid_atc_l3[c]) and c not in dups_q]
            if rel_p3:
                cand_p3 = np.array([i for i in atc_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)
                rel_p3_indices = [cid_to_idx[c] for c in rel_p3]
                rand_p3 = rand_hit_at_k(len(cand_p3), len(rel_p3_indices), 10)
                query_specs["P3_dn"].append((q_idx, cand_p3, rel_p3_indices, rand_p3))
                
    for p in protocols:
        print(f"Protocol {p}: {len(query_specs[p])} queries")
        
    # Real baseline values from on-disk CSVs
    def get_real_values(csv_path: Path, method: str):
        with open(csv_path, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["method"] == method:
                    if r["metric"] == "hit@10":
                        h = float(r["mean"])
                    elif r["metric"] == "mrr":
                        m = float(r["mean"])
                    elif r["metric"] == "enrichment":
                        e = float(r["mean"])
        return {"hit10": h, "mrr": m, "enrichment": e}

    real_data = {
        ("P1", "M3"): get_real_values(RESULTS_DIR / "metrics_C1_P1.csv", "M3"),
        ("P1", "M4"): get_real_values(RESULTS_DIR / "metrics_C1_P1.csv", "M4"),
        ("P2a_0.4", "M3"): get_real_values(RESULTS_DIR / "metrics_C1_P2a_0.4.csv", "M3"),
        ("P2a_0.4", "M4"): get_real_values(RESULTS_DIR / "metrics_C1_P2a_0.4.csv", "M4"),
        ("P3_dn", "M3"): get_real_values(RESULTS_DIR / "metrics_C1_P3_dn.csv", "M3"),
        ("P3_dn", "M4"): get_real_values(RESULTS_DIR / "metrics_C1_P3_dn.csv", "M4"),
    }

    # 2. Permutation Runs
    N_PERMS = 20
    
    # --- M3 Shuffle (Targets) ---
    print("\nRunning M3 Target Shuffle Control (20 permutations, seed 42)...")
    target_drugs_idx = [i for i, cid in enumerate(pool_cids) if len(drugs[cid]["targets"]) > 0]
    n_target_drugs = len(target_drugs_idx)
    orig_target_sets = [drugs[cid]["targets"] for cid in pool_cids]
    print(f"Drugs with targets: {n_target_drugs} / {N}")
    
    rng_t = np.random.RandomState(42)
    m3_perm_results = {p: {"hit10": [], "mrr": [], "enrichment": []} for p in protocols}
    
    for perm_idx in range(N_PERMS):
        perm_sub = rng_t.permutation(n_target_drugs)
        shuffled_targets = [set() for _ in range(N)]
        for i_sub, j_sub in enumerate(perm_sub):
            shuffled_targets[target_drugs_idx[i_sub]] = orig_target_sets[target_drugs_idx[j_sub]]
            
        T_perm = compute_jaccard_matrix(shuffled_targets, N)
        for p in protocols:
            res = evaluate_method_on_queries(T_perm, query_specs[p])
            for m in ["hit10", "mrr", "enrichment"]:
                m3_perm_results[p][m].append(res[m])
        print(f"  M3 Permutation {perm_idx + 1}/{N_PERMS} complete.")
        
    # --- M4 Shuffle (Side Effects) ---
    print("\nRunning M4 Side-Effect Shuffle Control (20 permutations, seed 42)...")
    se_drugs_idx = [i for i, cid in enumerate(pool_cids) if len(drugs[cid]["side_effects"]) > 0]
    n_se_drugs = len(se_drugs_idx)
    orig_se_sets = [drugs[cid]["side_effects"] for cid in pool_cids]
    print(f"Drugs with side effects: {n_se_drugs} / {N}")
    
    rng_se = np.random.RandomState(42)
    m4_perm_results = {p: {"hit10": [], "mrr": [], "enrichment": []} for p in protocols}
    
    for perm_idx in range(N_PERMS):
        perm_sub = rng_se.permutation(n_se_drugs)
        shuffled_se = [set() for _ in range(N)]
        for i_sub, j_sub in enumerate(perm_sub):
            shuffled_se[se_drugs_idx[i_sub]] = orig_se_sets[se_drugs_idx[j_sub]]
            
        E_perm = compute_jaccard_matrix(shuffled_se, N)
        for p in protocols:
            res = evaluate_method_on_queries(E_perm, query_specs[p])
            for m in ["hit10", "mrr", "enrichment"]:
                m4_perm_results[p][m].append(res[m])
        print(f"  M4 Permutation {perm_idx + 1}/{N_PERMS} complete.")
        
    # 3. Compile and Save Results
    out_csv = RESULTS_DIR / "shuffle_control_results.csv"
    rows = []
    
    for protocol in protocols:
        for method, perm_dict in [("M3", m3_perm_results), ("M4", m4_perm_results)]:
            real = real_data[(protocol, method)]
            for metric in ["hit10", "mrr", "enrichment"]:
                vals = perm_dict[protocol][metric]
                mean_shuf = float(np.mean(vals))
                min_shuf = float(np.min(vals))
                max_shuf = float(np.max(vals))
                real_val = float(real[metric])
                rows.append({
                    "protocol": protocol,
                    "method": method,
                    "metric": metric,
                    "real_value": f"{real_val:.4f}",
                    "shuffle_mean": f"{mean_shuf:.4f}",
                    "shuffle_min": f"{min_shuf:.4f}",
                    "shuffle_max": f"{max_shuf:.4f}",
                    "shuffle_range": f"[{min_shuf:.4f}, {max_shuf:.4f}]",
                    "n_permutations": N_PERMS,
                    "random_seed": 42
                })
                
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "protocol", "method", "metric", "real_value", "shuffle_mean",
            "shuffle_min", "shuffle_max", "shuffle_range", "n_permutations", "random_seed"
        ])
        writer.writeheader()
        writer.writerows(rows)
        
    print(f"\nShuffle control complete in {time.time() - t0:.2f}s.")
    print(f"Results saved to {out_csv}")

if __name__ == "__main__":
    main()
