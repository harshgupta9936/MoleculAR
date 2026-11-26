"""
experiments/step8_se_leakage.py

Implements Requirement N4:
(i) Side-effect leakage diagnostic:
    - Lists top 30 side-effect terms with highest information about relevance:
      * Ranking A: Raw co-occurrence count in relevant pairs (user explicit example in N4)
      * Ranking B: Mutual Information I(X_t; Y) with indication relevance
    - Flags terms matching disease or indication names (Exact gold indication match, Substring match, Clinical disease).
    - Saves to experiments/results/se_leakage_terms.csv.

(ii) Re-runs M4 and M6 with flagged terms removed on Protocol P1:
    - Original P1 baseline
    - Ablation 1: Top co-occurrence terms matching indications/diseases removed
    - Ablation 2: Top Mutual Information disease terms removed
    - Ablation 3: All exact gold-indication-matching terms removed
    - Saves metrics to experiments/results/se_leakage_ablation.csv.
"""

import csv
import sys
import time
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, precompute_matrices,
    compute_average_ranks, bootstrap_ci,
    DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

LOGS_DIR = ROOT / "experiments" / "logs"

def main():
    print("="*80)
    print("STEP 8: SIDE-EFFECT LEAKAGE CHECK (N4)")
    print("="*80)
    t0 = time.time()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load Primary Data
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    matrices = precompute_matrices(drugs, pool_cids)

    cid_to_idx = matrices["cid_to_idx"]
    S = matrices["S"]
    T = matrices["T"]
    E_orig = matrices["E"]
    penalties = matrices["penalties"]
    N = len(pool_cids)

    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}

    # Gold indication names
    gold_ind_names = set()
    with open(GOLD_PRIMARY_PATH, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if "indication_name" in r:
                gold_ind_names.add(r["indication_name"].strip().lower())

    cids_with_se = [c for c in cids_with_inds if len(drugs[c]["side_effects"]) > 0]
    n_se_drugs = len(cids_with_se)

    total_pairs_with_se = 0
    rel_pairs_with_se = 0
    term_cooccur_all = Counter()
    term_cooccur_rel = Counter()

    for i, q in enumerate(cids_with_se):
        q_inds = cid_ind_sets[q]
        q_ses = drugs[q]["side_effects"]
        for j in range(i + 1, n_se_drugs):
            c = cids_with_se[j]
            c_inds = cid_ind_sets[c]
            c_ses = drugs[c]["side_effects"]
            common = q_ses & c_ses
            is_rel = not q_inds.isdisjoint(c_inds)
            total_pairs_with_se += 1
            if is_rel:
                rel_pairs_with_se += 1
            for t in common:
                term_cooccur_all[t] += 1
                if is_rel:
                    term_cooccur_rel[t] += 1

    base_rate = rel_pairs_with_se / total_pairs_with_se

    term_stats = []
    for t, n_rel in term_cooccur_rel.items():
        n_all = term_cooccur_all[t]
        if n_all < 5:
            continue
        prec = n_rel / n_all
        lift = prec / base_rate

        p11 = n_rel / total_pairs_with_se
        p10 = (n_all - n_rel) / total_pairs_with_se
        p01 = (rel_pairs_with_se - n_rel) / total_pairs_with_se
        p00 = (total_pairs_with_se - rel_pairs_with_se - (n_all - n_rel)) / total_pairs_with_se
        px1 = n_all / total_pairs_with_se
        px0 = 1.0 - px1
        py1 = base_rate
        py0 = 1.0 - base_rate
        mi = 0.0
        for pxy, px, py in [(p11, px1, py1), (p10, px1, py0), (p01, px0, py1), (p00, px0, py0)]:
            if pxy > 0:
                mi += pxy * np.log2(pxy / (px * py))

        t_lower = t.lower()
        exact_match = t_lower in gold_ind_names
        sub_match = any(t_lower == ind or t_lower in ind for ind in gold_ind_names)

        # Clinical disease classification
        # Known clinical disease keywords
        disease_suffixes = ("itis", "osis", "emia", "aemia", "pathy", "oma", "syndrome", "disease", "infection", "disorder", "rupture", "colitis", "depression", "hypotension", "hypertension")
        is_disease_term = exact_match or sub_match or any(t_lower.endswith(suf) for suf in disease_suffixes) or ("pain" in t_lower) or ("shock" in t_lower)

        term_stats.append({
            "term": t,
            "n_rel": n_rel,
            "n_all": n_all,
            "prec": prec,
            "lift": lift,
            "mi": mi,
            "exact_indication_match": exact_match,
            "substring_indication_match": sub_match,
            "is_disease_or_indication": is_disease_term,
        })

    # Ranking 1: Raw co-occurrence in relevant pairs
    by_cooccur = sorted(term_stats, key=lambda x: x["n_rel"], reverse=True)[:30]
    # Ranking 2: Mutual Information
    by_mi = sorted(term_stats, key=lambda x: x["mi"], reverse=True)[:30]

    # Save candidate terms CSV
    terms_csv = RESULTS_DIR / "se_leakage_terms.csv"
    with open(terms_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ranking_type", "rank", "term", "n_cooccur_relevant", "n_cooccur_total", "precision", "lift", "mutual_info_bits", "exact_gold_indication_match", "substring_indication_match", "flagged_disease_or_indication"])
        for idx, r in enumerate(by_cooccur, 1):
            w.writerow(["raw_cooccurrence_rel", idx, r["term"], r["n_rel"], r["n_all"], f"{r['prec']:.4f}", f"{r['lift']:.2f}", f"{r['mi']:.6f}", r["exact_indication_match"], r["substring_indication_match"], r["is_disease_or_indication"]])
        for idx, r in enumerate(by_mi, 1):
            w.writerow(["mutual_information", idx, r["term"], r["n_rel"], r["n_all"], f"{r['prec']:.4f}", f"{r['lift']:.2f}", f"{r['mi']:.6f}", r["exact_indication_match"], r["substring_indication_match"], r["is_disease_or_indication"]])
    print(f"Saved {terms_csv}")

    # (ii) Re-run M4 and M6 on P1 with flagged terms removed
    # Define term exclusion sets:
    # Set 1: Top 30 raw co-occurrence terms flagged as disease/indication
    flagged_cooccur = {r["term"] for r in by_cooccur if r["is_disease_or_indication"]}
    # Set 2: Top 30 MI terms flagged as disease/indication
    flagged_mi = {r["term"] for r in by_mi if r["is_disease_or_indication"]}
    # Set 3: All exact gold-indication-matching terms in pool
    flagged_exact = {r["term"] for r in term_stats if r["exact_indication_match"]}

    ablation_sets = [
        ("original", set()),
        ("ablation_top30_cooccur_flagged", flagged_cooccur),
        ("ablation_top30_mi_flagged", flagged_mi),
        ("ablation_all_exact_indication_matches", flagged_exact),
    ]

    p1_queries = []
    for q_cid in cids_with_inds:
        q_idx = cid_to_idx[q_cid]
        q_inds = cid_ind_sets[q_cid]
        rel_P1 = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(cid_ind_sets[c])]
        if rel_P1:
            p1_queries.append({
                "q_cid": q_cid,
                "q_idx": q_idx,
                "cand_indices": np.array([i for i in range(N) if i != q_idx], dtype=np.int32),
                "rel_indices": np.array([cid_to_idx[c] for c in rel_P1], dtype=np.int32),
            })

    ablation_rows = []

    for set_name, terms_to_remove in ablation_sets:
        print(f"\n--- Running Ablation: {set_name} (Removing {len(terms_to_remove)} SE terms) ---")

        # Build ablated E matrix
        E_abl = np.zeros((N, N), dtype=np.float64)
        for i in range(N):
            se_i = drugs[pool_cids[i]]["side_effects"] - terms_to_remove
            if not se_i:
                continue
            for j in range(i + 1, N):
                se_j = drugs[pool_cids[j]]["side_effects"] - terms_to_remove
                if not se_j:
                    continue
                union_len = len(se_i | se_j)
                if union_len > 0:
                    val = len(se_i & se_j) / union_len
                    E_abl[i, j] = val
                    E_abl[j, i] = val

        # Evaluate M4 and M6 across all 1,661 P1 queries
        m4_hits, m4_mrrs = [], []
        m6_hits, m6_mrrs = [], []

        for q in p1_queries:
            q_idx = q["q_idx"]
            cands = q["cand_indices"]
            rels = q["rel_indices"]

            s_row = S[q_idx]
            t_row = T[q_idx]
            e_row = E_abl[q_idx]

            # M4: E alone
            m4_ranks = compute_average_ranks(e_row, cands)
            m4_rmin = min(m4_ranks[r] for r in rels)
            m4_hits.append(1.0 if m4_rmin <= 10.0 else 0.0)
            m4_mrrs.append(1.0 / m4_rmin)

            # M6: 0.4 S + 0.4 T + 0.2 E
            m6_scores = 0.4 * s_row + 0.4 * t_row + 0.2 * e_row
            m6_ranks = compute_average_ranks(m6_scores, cands)
            m6_rmin = min(m6_ranks[r] for r in rels)
            m6_hits.append(1.0 if m6_rmin <= 10.0 else 0.0)
            m6_mrrs.append(1.0 / m6_rmin)

        # Compute bootstrap CIs
        m4_h10_m, m4_h10_l, m4_h10_h = bootstrap_ci(np.array(m4_hits), n_boot=1000, seed=42)
        m4_mrr_m, m4_mrr_l, m4_mrr_h = bootstrap_ci(np.array(m4_mrrs), n_boot=1000, seed=42)

        m6_h10_m, m6_h10_l, m6_h10_h = bootstrap_ci(np.array(m6_hits), n_boot=1000, seed=42)
        m6_mrr_m, m6_mrr_l, m6_mrr_h = bootstrap_ci(np.array(m6_mrrs), n_boot=1000, seed=42)

        n_q = len(p1_queries)
        ablation_rows.append({"ablation_condition": set_name, "n_terms_removed": len(terms_to_remove), "method": "M4", "n_queries": n_q, "hit10": m4_h10_m, "hit10_ci_low": m4_h10_l, "hit10_ci_high": m4_h10_h, "mrr": m4_mrr_m, "mrr_ci_low": m4_mrr_l, "mrr_ci_high": m4_mrr_h})
        ablation_rows.append({"ablation_condition": set_name, "n_terms_removed": len(terms_to_remove), "method": "M6", "n_queries": n_q, "hit10": m6_h10_m, "hit10_ci_low": m6_h10_l, "hit10_ci_high": m6_h10_h, "mrr": m6_mrr_m, "mrr_ci_low": m6_mrr_l, "mrr_ci_high": m6_mrr_h})

        print(f"  M4: Hit@10 = {m4_h10_m:.4f} [{m4_h10_l:.4f}, {m4_h10_h:.4f}], MRR = {m4_mrr_m:.4f} [{m4_mrr_l:.4f}, {m4_mrr_h:.4f}]")
        print(f"  M6: Hit@10 = {m6_h10_m:.4f} [{m6_h10_l:.4f}, {m6_h10_h:.4f}], MRR = {m6_mrr_m:.4f} [{m6_mrr_l:.4f}, {m6_mrr_h:.4f}]")

    ablation_csv = RESULTS_DIR / "se_leakage_ablation.csv"
    with open(ablation_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ablation_condition", "n_terms_removed", "method", "n_queries", "hit10", "hit10_ci_low", "hit10_ci_high", "mrr", "mrr_ci_low", "mrr_ci_high"])
        w.writeheader()
        w.writerows(ablation_rows)
    print(f"\nSaved {ablation_csv} in {time.time() - t0:.2f}s")

if __name__ == "__main__":
    main()
