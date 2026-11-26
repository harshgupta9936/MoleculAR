"""
experiments/step19_known_cases_db_a.py

Requirement F3: Known-Case Check on DB-A
- Evaluates the 8 canonical benchmark cases on DB-A:
    1. Sildenafil -> Vardenafil
    2. Minoxidil -> Diazoxide
    3. Metformin -> Phenformin
    4. Thalidomide -> Lenalidomide
    5. Duloxetine -> Fluoxetine
    6. Gabapentin -> Pregabalin
    7. Bupropion -> Varenicline
    8. Raloxifene -> Fulvestrant
- Evaluates under Methods: M2 (Structural), M6 (Unpenalized Composite), M7 (Hubness-Penalized Composite), M8 (Popularity).
- Reports: Target drug rank, candidate pool size M, and top-10 list for M7.
- Strictly reports ranks and data without evaluative commentary.
- Outputs saved to experiments/results/known_cases_db_a.csv and printed.
"""

import sys
import csv
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, precompute_matrices,
    compute_average_ranks, DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

PAIRS = [
    ("Sildenafil", "CHEMBL1737", "Vardenafil", "CHEMBL1339"),
    ("Minoxidil", "CHEMBL802", "Diazoxide", "CHEMBL181"),
    ("Metformin", "CHEMBL1431", "Phenformin", "CHEMBL1528839"),
    ("Thalidomide", "CHEMBL468", "Lenalidomide", "CHEMBL6068337"),
    ("Duloxetine", "CHEMBL1175", "Fluoxetine", "CHEMBL1201082"),
    ("Gabapentin", "CHEMBL940", "Pregabalin", "CHEMBL1059"),
    ("Bupropion", "CHEMBL1201735", "Varenicline", "CHEMBL1076903"),
    ("Raloxifene", "CHEMBL1116", "Fulvestrant", "CHEMBL1358"),
]

def run_known_cases():
    print("=" * 85)
    print("STEP F3: KNOWN-CASE BENCHMARK CHECK ON DB-A")
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

    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    pop_scores = np.array([len(drug_inds[c]) for c in pool_cids], dtype=np.float64)

    # 1. Evaluate candidate pool c != q (M = 2,191)
    case_summary = []
    m7_top10_records = []

    for q_name, q_cid, t_name, t_cid in PAIRS:
        assert q_cid in cid_to_idx, f"Query CID {q_cid} ({q_name}) not in pool!"
        assert t_cid in cid_to_idx, f"Target CID {t_cid} ({t_name}) not in pool!"

        q_idx = cid_to_idx[q_cid]
        t_idx = cid_to_idx[t_cid]

        # Candidate pool: all pool drugs except query
        cand_indices = np.array([i for i in all_pool_indices if i != q_idx], dtype=np.int32)
        M_q = len(cand_indices)

        s = S[q_idx]
        t = T[q_idx]
        e = E[q_idx]
        m6 = 0.4 * s + 0.4 * t + 0.2 * e
        m7 = m6 / penalties
        m8 = pop_scores

        score_dict = {
            "M2": s,
            "M6": m6,
            "M7": m7,
            "M8": m8,
        }

        ranks = {}
        for m, sc in score_dict.items():
            r_map = compute_average_ranks(sc, cand_indices)
            ranks[m] = r_map[t_idx]

        case_summary.append({
            "query_drug": q_name,
            "query_chembl_id": q_cid,
            "target_drug": t_name,
            "target_chembl_id": t_cid,
            "candidate_pool_size": M_q,
            "rank_M2": ranks["M2"],
            "rank_M6": ranks["M6"],
            "rank_M7": ranks["M7"],
            "rank_M8": ranks["M8"],
        })

        # Top 10 for M7
        cand_m7_scores = [(pool_cids[i], drugs[pool_cids[i]]["name"], m7[i]) for i in cand_indices]
        # Sort descending by score, ascending by cid for tie-breaking in presentation
        cand_m7_scores.sort(key=lambda x: (-x[2], x[0]))

        r_map_m7 = compute_average_ranks(m7, cand_indices)
        for rank_pos, (c_cid, c_name, c_score) in enumerate(cand_m7_scores[:10], start=1):
            c_idx = cid_to_idx[c_cid]
            m7_top10_records.append({
                "query_drug": q_name,
                "query_chembl_id": q_cid,
                "position": rank_pos,
                "average_rank": r_map_m7[c_idx],
                "candidate_name": c_name,
                "candidate_chembl_id": c_cid,
                "m7_score": round(float(c_score), 6),
            })

    # Save summary table
    df_summary = pd.DataFrame(case_summary)
    summary_path = RESULTS_DIR / "known_cases_db_a.csv"
    df_summary.to_csv(summary_path, index=False)
    print(f"Saved known cases summary to: {summary_path}")

    # Save top 10 table
    df_top10 = pd.DataFrame(m7_top10_records)
    top10_path = RESULTS_DIR / "known_cases_m7_top10.csv"
    df_top10.to_csv(top10_path, index=False)
    print(f"Saved known cases M7 top 10 to: {top10_path}")

    print("\n--- Summary of 8 Known-Case Ranks on DB-A ---")
    print(df_summary.to_string())

if __name__ == "__main__":
    run_known_cases()
