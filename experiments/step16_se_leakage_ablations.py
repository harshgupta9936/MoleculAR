"""
experiments/step16_se_leakage_ablations.py

Requirement N12:
Leakage Ablations on Approved Term Lists for Protocols P1, P3-dn, and P4:
- Methods: M4 (Side-Effect Jaccard) and M6 (Chemical Structure + Side-Effect Fusion)
- Conditions:
  (A) Remove 7 indication-like terms:
      Intermittent claudication, Autoimmune disorder, Uterine leiomyoma,
      Acute leukaemia, Drug dependence, Withdrawal syndrome, Hyperprolactinaemia
  (B) Remove all 50 terms in experiments/results/se_top50_lift.csv
  (C) Remove (A) and (B) together (union of A and B)
- Reports:
  - Number of unique side-effect terms removed
  - Number of drugs losing all their side effects
  - Hit@10 and MRR with 95% bootstrap CIs (1,000 resamples, seed 42)
  - Comparison next to original values from existing CSVs
- Saves to experiments/results/se_leakage_approved_ablations.csv.
- Strictly no drawing conclusions.
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
    compute_average_ranks, bootstrap_ci,
    DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)
from experiments.step10_shuffle_control import compute_jaccard_matrix
from experiments.step11_p4_protocol import load_atc_level1

def main():
    print("=" * 80)
    print("STEP 16: APPROVED SIDE-EFFECT LEAKAGE ABLATIONS (REQUIREMENT N12)")
    print("=" * 80)
    t0 = time.time()
    
    # 1. Load Data
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    N = len(pool_cids)
    cid_to_idx = {cid: i for i, cid in enumerate(pool_cids)}
    
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    cid_atc_l3 = load_atc_l3(set(pool_cids))
    cid_atc_l1, _ = load_atc_level1(set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids)
    
    # Precompute S matrix
    fps = [drugs[cid]["fingerprint"] for cid in pool_cids]
    from rdkit.Chem import DataStructs
    S = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        if fps[i] is not None:
            sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
            S[i] = np.array(sims, dtype=np.float64)

    # Precompute T matrix
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
            
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    pool_with_atc_l3 = {cid for cid in pool_cids if len(cid_atc_l3[cid]) > 0}
    atc_l3_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc_l3], dtype=np.int32)
    pool_with_atc_l1 = {cid for cid in pool_cids if len(cid_atc_l1[cid]) > 0}
    atc_l1_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc_l1], dtype=np.int32)

    # Define protocols
    protocols = ["P1", "P3_dn", "P4"]
    proto_queries = {p: [] for p in protocols}

    for q_cid in cids_with_inds:
        q_idx = cid_to_idx[q_cid]
        q_inds = cid_ind_sets[q_cid]
        dups_q = duplicates_of.get(q_cid, set())
        
        rel_P1 = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(cid_ind_sets[c])]
        if not rel_P1:
            continue
            
        # P1
        cand_p1 = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)
        rel_p1_indices = np.array([cid_to_idx[c] for c in rel_P1], dtype=np.int32)
        proto_queries["P1"].append((q_idx, cand_p1, rel_p1_indices))
        
        # P3_dn
        if q_cid in pool_with_atc_l3:
            q_atc = cid_atc_l3[q_cid]
            rel_p3 = [c for c in rel_P1 if c in pool_with_atc_l3 and q_atc.isdisjoint(cid_atc_l3[c]) and c not in dups_q]
            if rel_p3:
                cand_p3 = np.array([i for i in atc_l3_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)
                rel_p3_indices = np.array([cid_to_idx[c] for c in rel_p3], dtype=np.int32)
                proto_queries["P3_dn"].append((q_idx, cand_p3, rel_p3_indices))
                
        # P4
        if q_cid in pool_with_atc_l1:
            q_l1 = cid_atc_l1[q_cid]
            rel_p4 = [c for c in rel_P1 if c in pool_with_atc_l1 and q_l1.isdisjoint(cid_atc_l1[c]) and c not in dups_q]
            if rel_p4:
                cand_p4 = np.array([i for i in atc_l1_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)
                rel_p4_indices = np.array([cid_to_idx[c] for c in rel_p4], dtype=np.int32)
                proto_queries["P4"].append((q_idx, cand_p4, rel_p4_indices))

    # Define term lists
    list_a_terms = {
        "Intermittent claudication",
        "Autoimmune disorder",
        "Uterine leiomyoma",
        "Acute leukaemia",
        "Drug dependence",
        "Withdrawal syndrome",
        "Hyperprolactinaemia"
    }
    
    list_b_terms = set()
    with open(RESULTS_DIR / "se_top50_lift.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            list_b_terms.add(r["term"].strip())
            
    list_c_terms = list_a_terms | list_b_terms

    conditions = [
        ("Original (None removed)", set()),
        ("(A) 7 Indication-like terms removed", list_a_terms),
        ("(B) 50 Top Lift terms removed", list_b_terms),
        ("(C) (A) + (B) Together removed", list_c_terms)
    ]

    # Baseline original values from existing CSVs
    def get_orig(csv_name, method, metric):
        with open(RESULTS_DIR / csv_name, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["method"] == method and r["metric"] == metric:
                    return float(r["mean"]), float(r["ci_low"]), float(r["ci_high"])
        return None, None, None

    orig_baselines = {
        ("P1", "M4", "hit@10"): get_orig("metrics_C1_P1.csv", "M4", "hit@10"),
        ("P1", "M4", "mrr"): get_orig("metrics_C1_P1.csv", "M4", "mrr"),
        ("P1", "M6", "hit@10"): get_orig("metrics_C1_P1.csv", "M6", "hit@10"),
        ("P1", "M6", "mrr"): get_orig("metrics_C1_P1.csv", "M6", "mrr"),
        ("P3_dn", "M4", "hit@10"): get_orig("metrics_C1_P3_dn.csv", "M4", "hit@10"),
        ("P3_dn", "M4", "mrr"): get_orig("metrics_C1_P3_dn.csv", "M4", "mrr"),
        ("P3_dn", "M6", "hit@10"): get_orig("metrics_C1_P3_dn.csv", "M6", "hit@10"),
        ("P3_dn", "M6", "mrr"): get_orig("metrics_C1_P3_dn.csv", "M6", "mrr"),
        ("P4", "M4", "hit@10"): get_orig("metrics_C1_P4.csv", "M4", "hit@10"),
        ("P4", "M4", "mrr"): get_orig("metrics_C1_P4.csv", "M4", "mrr"),
        ("P4", "M6", "hit@10"): get_orig("metrics_C1_P4.csv", "M6", "hit@10"),
        ("P4", "M6", "mrr"): get_orig("metrics_C1_P4.csv", "M6", "mrr"),
    }

    out_rows = []
    
    orig_drugs_with_se = sum(1 for cid in pool_cids if len(drugs[cid]["side_effects"]) > 0)
    print(f"Original pool drugs with side effects: {orig_drugs_with_se}")

    for cond_name, remove_set in conditions:
        n_terms_removed = len(remove_set)
        # Create filtered side effect sets
        filtered_ses = []
        drugs_lost_all_se = 0
        for cid in pool_cids:
            orig_set = drugs[cid]["side_effects"]
            new_set = {t for t in orig_set if t not in remove_set}
            if len(orig_set) > 0 and len(new_set) == 0:
                drugs_lost_all_se += 1
            filtered_ses.append(new_set)
            
        print(f"\nCondition: {cond_name}")
        print(f"  Terms removed: {n_terms_removed}")
        print(f"  Drugs losing all side effects: {drugs_lost_all_se}")
        
        # Compute ablated E matrix
        E_abl = compute_jaccard_matrix(filtered_ses, N)
        
        # Evaluate for P1, P3_dn, P4
        for p in protocols:
            queries = proto_queries[p]
            n_queries = len(queries)
            
            # M4 scores
            m4_hits = []
            m4_mrrs = []
            # M6 scores
            m6_hits = []
            m6_mrrs = []
            
            for q_idx, cands, rels in queries:
                e_scores = E_abl[q_idx]
                s_scores = S[q_idx]
                t_scores = T[q_idx]
                m6_scores = 0.4 * s_scores + 0.4 * t_scores + 0.2 * e_scores
                
                # M4 ranks
                ranks_m4 = compute_average_ranks(e_scores, cands)
                rel_ranks_m4 = [ranks_m4[c] for c in rels]
                r_min_m4 = min(rel_ranks_m4)
                m4_hits.append(1.0 if r_min_m4 <= 10.0 else 0.0)
                m4_mrrs.append(1.0 / r_min_m4)
                
                # M6 ranks
                ranks_m6 = compute_average_ranks(m6_scores, cands)
                rel_ranks_m6 = [ranks_m6[c] for c in rels]
                r_min_m6 = min(rel_ranks_m6)
                m6_hits.append(1.0 if r_min_m6 <= 10.0 else 0.0)
                m6_mrrs.append(1.0 / r_min_m6)
                
            # Bootstrap CIs (1,000 resamples, seed 42)
            h10_m4_m, h10_m4_l, h10_m4_u = bootstrap_ci(np.array(m4_hits), n_boot=1000, seed=42)
            mrr_m4_m, mrr_m4_l, mrr_m4_u = bootstrap_ci(np.array(m4_mrrs), n_boot=1000, seed=42)
            
            h10_m6_m, h10_m6_l, h10_m6_u = bootstrap_ci(np.array(m6_hits), n_boot=1000, seed=42)
            mrr_m6_m, mrr_m6_l, mrr_m6_u = bootstrap_ci(np.array(m6_mrrs), n_boot=1000, seed=42)
            
            for m, (h_m, h_l, h_u), (m_m, m_l, m_u) in [
                ("M4", (h10_m4_m, h10_m4_l, h10_m4_u), (mrr_m4_m, mrr_m4_l, mrr_m4_u)),
                ("M6", (h10_m6_m, h10_m6_l, h10_m6_u), (mrr_m6_m, mrr_m6_l, mrr_m6_u))
            ]:
                orig_h_m, orig_h_l, orig_h_u = orig_baselines[(p, m, "hit@10")]
                orig_m_m, orig_m_l, orig_m_u = orig_baselines[(p, m, "mrr")]
                
                out_rows.append({
                    "protocol": p,
                    "condition": cond_name,
                    "n_terms_removed": n_terms_removed,
                    "drugs_lost_all_se": drugs_lost_all_se,
                    "method": m,
                    "hit10_mean": f"{h_m:.4f}",
                    "hit10_ci": f"{h_m:.4f} [{h_l:.4f}, {h_u:.4f}]",
                    "orig_hit10_ci": f"{orig_h_m:.4f} [{orig_h_l:.4f}, {orig_h_u:.4f}]",
                    "mrr_mean": f"{m_m:.4f}",
                    "mrr_ci": f"{m_m:.4f} [{m_l:.4f}, {m_u:.4f}]",
                    "orig_mrr_ci": f"{orig_m_m:.4f} [{orig_m_l:.4f}, {orig_m_u:.4f}]",
                    "n_queries": n_queries,
                })

    out_csv = RESULTS_DIR / "se_leakage_approved_ablations.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "protocol", "condition", "n_terms_removed", "drugs_lost_all_se", "method",
            "hit10_mean", "hit10_ci", "orig_hit10_ci",
            "mrr_mean", "mrr_ci", "orig_mrr_ci", "n_queries"
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)
        
    print(f"\nAblations complete in {time.time() - t0:.2f}s.")
    print(f"Results saved to {out_csv}")

if __name__ == "__main__":
    main()
