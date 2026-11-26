"""
experiments/step12_se_lift_terms.py

Requirement N9:
Produce the top 50 side-effect terms ranked by lift with minimum support (at least 20 drugs with the term),
with columns:
- term
- n drugs
- n relevant pair co-occurrences
- lift

Saves to: experiments/results/se_top50_lift.csv
"""

import sys
import csv
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs,
    DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

def main():
    print("=" * 80)
    print("STEP 12: TOP 50 SIDE-EFFECT TERMS RANKED BY LIFT (REQUIREMENT N9)")
    print("=" * 80)
    
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Load Data
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    
    # Drugs in evaluation set: have >= 1 indication AND >= 1 side effect
    cids_with_se = [c for c in cids_with_inds if len(drugs[c]["side_effects"]) > 0]
    n_eval_drugs = len(cids_with_se)
    print(f"Total pool drugs: {len(pool_cids)}")
    print(f"Drugs with indications: {len(cids_with_inds)}")
    print(f"Drugs with indications and side effects: {n_eval_drugs}")
    
    # Drug count for each term across all pool drugs, and across eval drugs
    term_drug_count_pool = Counter()
    for cid in pool_cids:
        for t in drugs[cid]["side_effects"]:
            term_drug_count_pool[t] += 1
            
    term_drug_count_eval = Counter()
    for cid in cids_with_se:
        for t in drugs[cid]["side_effects"]:
            term_drug_count_eval[t] += 1

    # Pairwise co-occurrences in evaluation set
    total_pairs = 0
    rel_pairs = 0
    term_cooccur_all = Counter()
    term_cooccur_rel = Counter()
    
    for i in range(n_eval_drugs):
        q = cids_with_se[i]
        q_inds = cid_ind_sets[q]
        q_ses = drugs[q]["side_effects"]
        for j in range(i + 1, n_eval_drugs):
            c = cids_with_se[j]
            c_inds = cid_ind_sets[c]
            c_ses = drugs[c]["side_effects"]
            common = q_ses & c_ses
            is_rel = not q_inds.isdisjoint(c_inds)
            
            total_pairs += 1
            if is_rel:
                rel_pairs += 1
            for t in common:
                term_cooccur_all[t] += 1
                if is_rel:
                    term_cooccur_rel[t] += 1
                    
    base_rate = rel_pairs / total_pairs
    print(f"Total drug pairs: {total_pairs}, Relevant drug pairs: {rel_pairs}, Base rate: {base_rate:.6f}")
    
    # Filter terms with minimum support: at least 20 drugs with the term
    # We use term_drug_count_eval >= 20 (and verify pool count >= 20)
    term_rows = []
    for t, n_drugs in term_drug_count_eval.items():
        if n_drugs < 20:
            continue
        n_all_cooccur = term_cooccur_all[t]
        n_rel_cooccur = term_cooccur_rel[t]
        if n_all_cooccur == 0:
            continue
        precision = n_rel_cooccur / n_all_cooccur
        lift = precision / base_rate
        term_rows.append({
            "term": t,
            "n_drugs": n_drugs,
            "n_drugs_pool": term_drug_count_pool[t],
            "n_relevant_pair_cooccurrences": n_rel_cooccur,
            "n_total_pair_cooccurrences": n_all_cooccur,
            "lift": lift
        })
        
    # Sort by lift descending
    term_rows.sort(key=lambda x: x["lift"], reverse=True)
    top50 = term_rows[:50]
    
    out_csv = RESULTS_DIR / "se_top50_lift.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "rank", "term", "n_drugs", "n_relevant_pair_cooccurrences", "lift"
        ])
        writer.writeheader()
        for rank, r in enumerate(top50, 1):
            writer.writerow({
                "rank": rank,
                "term": r["term"],
                "n_drugs": r["n_drugs"],
                "n_relevant_pair_cooccurrences": r["n_relevant_pair_cooccurrences"],
                "lift": f"{r['lift']:.4f}"
            })
            
    print(f"Top 50 terms written to {out_csv}")
    for r in top50[:10]:
        print(f"  {r['term']}: n_drugs={r['n_drugs']}, n_rel={r['n_relevant_pair_cooccurrences']}, lift={r['lift']:.4f}")

if __name__ == "__main__":
    main()
