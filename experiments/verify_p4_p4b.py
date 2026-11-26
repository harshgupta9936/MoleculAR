"""
experiments/verify_p4_p4b.py

Step F1: Verification for Protocols P4 and P4b
- Brute-force rank calculation for M2 to M8 on:
    * 10 random qualifying queries (seed 42)
    * 5 queries with no targets (seed 42)
    * 5 queries with no side effects (seed 42)
- Evaluated on both Protocol P4 and Protocol P4b.
- Exact tie rule written directly in the verifier with NO shared helper:
    rank = 1.0 + N_greater + 0.5 * N_equal
- Asserts brute-force ranks match engine ranks within 1e-12.
- Pytest function test_step14_p1_reproduction() verifying that step14 evaluation logic
  under P1-equivalent settings reproduces experiments/results/metrics_C1_P1.csv within 1e-9.
"""

import sys
import random
import sqlite3
import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_duplicates, precompute_matrices,
    compute_average_ranks, DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)
from experiments.step11_p4_protocol import load_atc_level1

def run_verification(verbose=True):
    log_lines = []
    def log(msg=""):
        if verbose:
            print(msg)
        log_lines.append(str(msg))

    log("=" * 85)
    log("STEP F1: INDEPENDENT BRUTE-FORCE VERIFICATION FOR PROTOCOLS P4 AND P4b")
    log("=" * 85)

    # 1. Load Main Engine matrices and data
    main_drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(main_drugs.keys()))
    main_matrices = precompute_matrices(main_drugs, pool_cids)
    cid_to_idx = main_matrices["cid_to_idx"]
    N = len(pool_cids)

    cid_atc_l1, _ = load_atc_level1(set(pool_cids))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids)

    pop_scores = np.array([len(drug_inds[c]) for c in pool_cids], dtype=np.float64)

    # 2. Independent Raw Data Loading directly from SQLite (zero shared code)
    conn = sqlite3.connect(DB_A_PATH.resolve().as_uri() + "?mode=ro", uri=True)
    cur = conn.cursor()

    cur.execute("SELECT chembl_id, fingerprint_blob, hubness_percentile FROM drugs WHERE trivial_category IS NULL")
    raw_drugs = {}
    for cid, fp_blob, hub in cur.fetchall():
        if cid in pool_cids:
            bits = set()
            if fp_blob:
                for idx, ch in enumerate(fp_blob):
                    if ch == "1":
                        bits.add(idx)
            raw_drugs[cid] = {
                "bits": bits,
                "hub": hub if hub is not None else 0.0,
                "targets": set(),
                "side_effects": set(),
                "indications": set(drug_inds.get(cid, set())),
            }

    cur.execute("SELECT chembl_id, uniprot_id FROM targets")
    for cid, uid in cur.fetchall():
        if cid in raw_drugs and uid:
            raw_drugs[cid]["targets"].add(uid)

    cur.execute("SELECT chembl_id, side_effect_name FROM side_effects")
    for cid, se in cur.fetchall():
        if cid in raw_drugs and se:
            raw_drugs[cid]["side_effects"].add(se)
    conn.close()

    # 3. Identify qualifying queries for P4 / P4b
    qualifying_queries = []
    for q in pool_cids:
        if q not in cid_atc_l1 or not cid_atc_l1[q] or not drug_inds[q]:
            continue
        q_l1 = cid_atc_l1[q]
        q_inds = drug_inds[q]
        dups_q = duplicates_of.get(q, set())
        rels = [
            c for c in pool_cids
            if c != q and c not in dups_q and c in cid_atc_l1 and cid_atc_l1[c]
            and not q_inds.isdisjoint(drug_inds[c]) and q_l1.isdisjoint(cid_atc_l1[c])
        ]
        if rels:
            qualifying_queries.append(q)

    log(f"Total Qualifying Queries for P4/P4b: N = {len(qualifying_queries)}")

    # 4. Sample Query Groups (seed 42)
    # (a) 10 random queries
    rng_rand = random.Random(42)
    sampled_random = sorted(rng_rand.sample(qualifying_queries, 10))

    # (b) 5 queries with no targets
    empty_t_pool = sorted([q for q in qualifying_queries if len(raw_drugs[q]["targets"]) == 0])
    rng_t = random.Random(42)
    sampled_empty_t = sorted(rng_t.sample(empty_t_pool, 5))

    # (c) 5 queries with no side effects
    empty_se_pool = sorted([q for q in qualifying_queries if len(raw_drugs[q]["side_effects"]) == 0])
    rng_se = random.Random(42)
    sampled_empty_se = sorted(rng_se.sample(empty_se_pool, 5))

    query_groups = [
        ("Random Qualifying Queries (n=10)", sampled_random),
        ("Queries with No Targets (n=5)", sampled_empty_t),
        ("Queries with No Side Effects (n=5)", sampled_empty_se),
    ]

    atc_pool = sorted([c for c in pool_cids if c in cid_atc_l1 and len(cid_atc_l1[c]) > 0])
    methods = ["M2", "M3", "M4", "M5", "M6", "M7", "M8"]

    for protocol_name in ["P4", "P4b"]:
        log("\n" + "=" * 80)
        log(f"VERIFYING PROTOCOL: {protocol_name}")
        log("=" * 80)

        for group_name, query_list in query_groups:
            log(f"\n--- Group: {group_name} ---")
            for q_cid in query_list:
                q_raw = raw_drugs[q_cid]
                q_idx = cid_to_idx[q_cid]
                q_l1 = cid_atc_l1[q_cid]
                dups_q = duplicates_of.get(q_cid, set())

                # Candidate pool definition
                if protocol_name == "P4":
                    cand_cids = [c for c in atc_pool if c != q_cid and c not in dups_q]
                else: # P4b: remove sharing ATC-L1
                    cand_cids = [
                        c for c in atc_pool
                        if c != q_cid and c not in dups_q and q_l1.isdisjoint(cid_atc_l1[c])
                    ]

                cand_indices = np.array([cid_to_idx[c] for c in cand_cids], dtype=np.int32)
                log(f"Query {q_cid} ({main_drugs[q_cid]['name']}): candidates={len(cand_cids)}, targets={len(q_raw['targets'])}, SEs={len(q_raw['side_effects'])}")

                # 1. Independent Brute Force Scores
                bf_scores = {m: {} for m in methods}
                for c_cid in cand_cids:
                    c_raw = raw_drugs[c_cid]

                    # S: bitwise Tanimoto
                    inter_bits = len(q_raw["bits"] & c_raw["bits"])
                    union_bits = len(q_raw["bits"] | c_raw["bits"])
                    s_sim = inter_bits / union_bits if union_bits > 0 else 0.0

                    # T: UniProt set Jaccard
                    inter_t = len(q_raw["targets"] & c_raw["targets"])
                    union_t = len(q_raw["targets"] | c_raw["targets"])
                    t_sim = inter_t / union_t if union_t > 0 else 0.0

                    # E: SE set Jaccard
                    inter_se = len(q_raw["side_effects"] & c_raw["side_effects"])
                    union_se = len(q_raw["side_effects"] | c_raw["side_effects"])
                    se_sim = inter_se / union_se if union_se > 0 else 0.0

                    # Hub penalty
                    hub_pct = c_raw["hub"]
                    hub_pen = 1.0 + (hub_pct / 100.0) ** 2

                    bf_scores["M2"][c_cid] = s_sim
                    bf_scores["M3"][c_cid] = t_sim
                    bf_scores["M4"][c_cid] = se_sim
                    bf_scores["M5"][c_cid] = 0.5 * s_sim + 0.5 * t_sim
                    m6_raw = 0.4 * s_sim + 0.4 * t_sim + 0.2 * se_sim
                    bf_scores["M6"][c_cid] = m6_raw
                    bf_scores["M7"][c_cid] = m6_raw / hub_pen
                    bf_scores["M8"][c_cid] = float(len(c_raw["indications"]))

                # 2. Independent Brute Force Ranks (Tie rule written directly with NO helper)
                bf_ranks = {m: {} for m in methods}
                for m in methods:
                    m_scores = bf_scores[m]
                    for c_cid in cand_cids:
                        sc = m_scores[c_cid]
                        # Direct formula: 1.0 + N_greater + 0.5 * N_equal
                        n_gt = sum(1 for o_cid in cand_cids if m_scores[o_cid] > sc)
                        n_eq = sum(1 for o_cid in cand_cids if m_scores[o_cid] == sc) - 1
                        bf_ranks[m][c_cid] = 1.0 + n_gt + 0.5 * n_eq

                # 3. Engine Ranks via canonical formula
                s_row = main_matrices["S"][q_idx]
                t_row = main_matrices["T"][q_idx]
                e_row = main_matrices["E"][q_idx]
                penalties = main_matrices["penalties"]
                m6_engine = 0.4 * s_row + 0.4 * t_row + 0.2 * e_row

                engine_score_dict = {
                    "M2": s_row,
                    "M3": t_row,
                    "M4": e_row,
                    "M5": 0.5 * s_row + 0.5 * t_row,
                    "M6": m6_engine,
                    "M7": m6_engine / penalties,
                    "M8": pop_scores,
                }

                for m in methods:
                    eng_scores = engine_score_dict[m]
                    eng_ranks = compute_average_ranks(eng_scores, cand_indices)
                    for c_cid in cand_cids:
                        c_idx = cid_to_idx[c_cid]
                        bf_r = bf_ranks[m][c_cid]
                        eng_r = eng_ranks[c_idx]
                        diff = abs(bf_r - eng_r)
                        assert diff < 1e-12, f"Rank mismatch in {protocol_name} {m} for query {q_cid}, candidate {c_cid}: bf={bf_r}, engine={eng_r}, diff={diff}"
                log(f"  -> All methods M2-M8 matched engine ranks within 1e-12 ({len(cand_cids)} candidates verified).")

    log("\n" + "=" * 85)
    log("ALL 40 VERIFICATION QUERIES (20 for P4, 20 for P4b) PASSED EXACTLY.")
    log("=" * 85)

    with open(ROOT / "experiments" / "logs" / "verify_p4_p4b.log", "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines) + "\n")

    return True

def test_verify_p4_p4b_brute_force():
    assert run_verification(verbose=False)

def test_step14_p1_reproduction():
    """
    Checks that the step14 evaluation logic for P1-equivalent settings reproduces
    experiments/results/metrics_C1_P1.csv within 1e-9 (after rounding to file precision 6 dec).
    """
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    matrices = precompute_matrices(drugs, pool_cids)
    cid_to_idx = matrices["cid_to_idx"]
    gold_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))

    queries = [c for c in pool_cids if len(gold_inds[c]) > 0]
    all_pool_indices = np.arange(len(pool_cids), dtype=np.int32)
    pop_scores = np.array([len(gold_inds[c]) for c in pool_cids], dtype=np.float64)

    df_csv = pd.read_csv(RESULTS_DIR / "metrics_C1_P1.csv")

    res = {m: {"hit@10": [], "hit@20": [], "precision@10": [], "mrr": [], "percentile_rank": []}
           for m in ["M2", "M3", "M4", "M5", "M6", "M7", "M8"]}

    for q_cid in queries:
        q_idx = cid_to_idx[q_cid]
        q_inds = gold_inds[q_cid]
        rels = [c for c in queries if c != q_cid and not q_inds.isdisjoint(gold_inds[c])]
        if not rels:
            continue
        cands = np.array([i for i in all_pool_indices if i != q_idx], dtype=np.int32)
        rel_indices = np.array([cid_to_idx[c] for c in rels], dtype=np.int32)

        s = matrices["S"][q_idx]
        t = matrices["T"][q_idx]
        e = matrices["E"][q_idx]
        m6 = 0.4 * s + 0.4 * t + 0.2 * e
        scores = {
            "M2": s,
            "M3": t,
            "M4": e,
            "M5": 0.5 * s + 0.5 * t,
            "M6": m6,
            "M7": m6 / matrices["penalties"],
            "M8": pop_scores,
        }
        M_q = len(cands)
        for m in scores:
            ranks = compute_average_ranks(scores[m], cands)
            rel_ranks = [ranks[c] for c in rel_indices]
            min_r = min(rel_ranks)
            n_rel_top10 = sum(1 for r in rel_ranks if r <= 10.0)

            res[m]["hit@10"].append(1.0 if min_r <= 10.0 else 0.0)
            res[m]["hit@20"].append(1.0 if min_r <= 20.0 else 0.0)
            res[m]["precision@10"].append(n_rel_top10 / 10.0)
            res[m]["mrr"].append(1.0 / min_r)
            res[m]["percentile_rank"].append(min_r / M_q)

    # Assert reproduction within 1e-9 (rounded to 6 decimal places as stored in CSV)
    for m in res:
        for met in ["hit@10", "hit@20", "precision@10", "mrr", "percentile_rank"]:
            c_mean = np.mean(res[m][met])
            csv_val = df_csv[(df_csv["method"] == m) & (df_csv["metric"] == met)]["mean"].values[0]
            # Rounded to 6 decimals
            diff_round = abs(round(c_mean, 6) - csv_val)
            assert diff_round < 1e-9, f"Reproduction mismatch for {m} {met}: calculated={c_mean:.8f}, rounded={round(c_mean, 6):.6f}, csv={csv_val:.6f}, diff={diff_round}"

if __name__ == "__main__":
    run_verification(verbose=True)
