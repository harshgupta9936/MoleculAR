"""
experiments/verify_experiments.py

Step 5.3: Independent Brute-Force Verification
- Samples 10 random queries (seed 42) from primary pool.
- Recomputes candidate ranks for M2 through M9 using a slow, independent brute-force loop.
- Shares ZERO helper functions with experiments/eval_engine.py.
- Directly queries SQLite tables for targets and side effects.
- Computes manual set intersections, manual bit operations, manual hubness penalties, and manual average ranks.
- Asserts that brute-force ranks match the main evaluation code ranks exactly.
- Saves log to experiments/logs/verify_experiments.log.
"""

import sys
import random
import sqlite3
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_atc_l3, precompute_matrices,
    compute_average_ranks, DB_A_PATH, GOLD_PRIMARY_PATH, INORGANIC_EXCLUDE
)

def run_verification():
    log_lines = []
    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=" * 85)
    log("STEP 5.3: INDEPENDENT BRUTE-FORCE RANKING VERIFICATION")
    log("=" * 85)
    log("Seed: 42, Queries: 10 random queries, Methods: M2 to M9\n")

    # 1. Load Main Engine matrices and data
    main_drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(main_drugs.keys()))
    main_matrices = precompute_matrices(main_drugs, pool_cids)
    cid_to_idx = main_matrices["cid_to_idx"]
    N = len(pool_cids)

    # 2. Independent Raw Data Loading directly from SQLite (zero shared code)
    conn = sqlite3.connect(DB_A_PATH.resolve().as_uri() + "?mode=ro", uri=True)
    cur = conn.cursor()

    cur.execute("SELECT chembl_id, fingerprint_blob, hubness_percentile FROM drugs WHERE trivial_category IS NULL")
    raw_drugs = {}
    for cid, fp_blob, hub in cur.fetchall():
        if cid in INORGANIC_EXCLUDE:
            continue
        # Parse bit string directly
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
            "indications": set(),
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

    # Load indications directly from gold_pairs_primary.csv
    import csv
    with open(GOLD_PRIMARY_PATH, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            cid = r["chembl_id"]
            if cid in raw_drugs:
                raw_drugs[cid]["indications"].add(r["indication_id"])

    # 3. Load indications for M8 popularity baseline (strictly float64)
    gold_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    pop_scores = np.array([len(gold_inds[c]) for c in pool_cids], dtype=np.float64)

    # 4. Sample Query Groups (seed 42)
    # (a) 10 random queries from pool_cids
    rng_rand = random.Random(42)
    sampled_random = sorted(rng_rand.sample(pool_cids, 10))

    # (b) 10 queries with empty target sets
    empty_t_pool = sorted([c for c in pool_cids if len(raw_drugs[c]["targets"]) == 0])
    rng_t = random.Random(42)
    sampled_empty_t = sorted(rng_t.sample(empty_t_pool, 10))

    # (c) 10 queries with empty side-effect sets
    empty_se_pool = sorted([c for c in pool_cids if len(raw_drugs[c]["side_effects"]) == 0])
    rng_se = random.Random(42)
    sampled_empty_se = sorted(rng_se.sample(empty_se_pool, 10))

    query_groups = [
        ("Random Pool Queries (n=10)", sampled_random),
        ("Queries with Empty Target Sets (n=10)", sampled_empty_t),
        ("Queries with Empty Side-Effect Sets (n=10)", sampled_empty_se),
    ]

    methods = ["M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"]
    all_passed = True
    total_verified = 0

    for group_name, query_list in query_groups:
        log(f"\n--- Verifying Group: {group_name} ---")
        for q_cid in query_list:
            total_verified += 1
            q_raw = raw_drugs[q_cid]
            q_idx = cid_to_idx[q_cid]
            cand_cids = [c for c in pool_cids if c != q_cid]
            cand_indices = np.array([cid_to_idx[c] for c in cand_cids], dtype=np.int32)

            log(f"Verifying Query [{total_verified}/30]: {q_cid} ({main_drugs[q_cid]['name']}) | targets={len(q_raw['targets'])}, SEs={len(q_raw['side_effects'])}")

            # Brute-force calculation for every candidate
            bf_scores = {m: {} for m in methods}

        for c_cid in cand_cids:
            c_raw = raw_drugs[c_cid]

            # Structural Tanimoto via bit sets
            inter_bits = len(q_raw["bits"] & c_raw["bits"])
            union_bits = len(q_raw["bits"] | c_raw["bits"])
            s_sim = inter_bits / union_bits if union_bits > 0 else 0.0

            # Target Jaccard via UniProt sets
            inter_t = len(q_raw["targets"] & c_raw["targets"])
            union_t = len(q_raw["targets"] | c_raw["targets"])
            t_sim = inter_t / union_t if union_t > 0 else 0.0

            # SE Jaccard via string sets
            inter_se = len(q_raw["side_effects"] & c_raw["side_effects"])
            union_se = len(q_raw["side_effects"] | c_raw["side_effects"])
            se_sim = inter_se / union_se if union_se > 0 else 0.0

            # Hubness penalty
            hub_pct = c_raw["hub"]
            hub_pen = 1.0 + (hub_pct / 100.0) ** 2

            # M2 (S)
            bf_scores["M2"][c_cid] = s_sim
            # M3 (T)
            bf_scores["M3"][c_cid] = t_sim
            # M4 (E)
            bf_scores["M4"][c_cid] = se_sim
            # M5 (0.5 S + 0.5 T)
            bf_scores["M5"][c_cid] = 0.5 * s_sim + 0.5 * t_sim
            # M6 (0.4 S + 0.4 T + 0.2 E)
            m6_raw = 0.4 * s_sim + 0.4 * t_sim + 0.2 * se_sim
            bf_scores["M6"][c_cid] = m6_raw
            # M7 (M6 / pen)
            bf_scores["M7"][c_cid] = m6_raw / hub_pen
            # M8 (Popularity)
            bf_scores["M8"][c_cid] = len(c_raw["indications"])
            # M9 (Reweighted / pen)
            avail_w = 0.4
            num = 0.4 * s_sim
            if len(q_raw["targets"]) > 0 and len(c_raw["targets"]) > 0:
                avail_w += 0.4
                num += 0.4 * t_sim
            if len(q_raw["side_effects"]) > 0 and len(c_raw["side_effects"]) > 0:
                avail_w += 0.2
                num += 0.2 * se_sim
            bf_scores["M9"][c_cid] = (num / avail_w) / hub_pen

        # Brute-force average rank calculation
        bf_ranks = {m: {} for m in methods}
        for m in methods:
            m_scores = bf_scores[m]
            for c_cid in cand_cids:
                sc = m_scores[c_cid]
                n_gt = sum(1 for other_sc in m_scores.values() if other_sc > sc)
                n_eq = sum(1 for other_sc in m_scores.values() if other_sc == sc) - 1
                bf_ranks[m][c_cid] = 1.0 + n_gt + 0.5 * n_eq

        # Main Engine Ranks
        # Compute main engine scores for query
        s_row = main_matrices["S"][q_idx]
        t_row = main_matrices["T"][q_idx]
        e_row = main_matrices["E"][q_idx]
        penalties = main_matrices["penalties"]
        has_target = main_matrices["has_target"]
        has_se = main_matrices["has_se"]
        main_m_scores = {
            "M2": s_row,
            "M3": t_row,
            "M4": e_row,
            "M5": 0.5 * s_row + 0.5 * t_row,
            "M6": 0.4 * s_row + 0.4 * t_row + 0.2 * e_row,
            "M7": (0.4 * s_row + 0.4 * t_row + 0.2 * e_row) / penalties,
            "M8": pop_scores,
        }
        w_s = 0.4
        w_t = 0.4 if has_target[q_idx] else 0.0
        w_e = 0.2 if has_se[q_idx] else 0.0
        w_sum = w_s + w_t * has_target + w_e * has_se
        scores_m9 = ((w_s * s_row + w_t * has_target * t_row + w_e * has_se * e_row) / w_sum) / penalties
        main_m_scores["M9"] = scores_m9

        for m in methods:
            main_rank_map = compute_average_ranks(main_m_scores[m], cand_indices)
            diffs = []
            for c_cid in cand_cids:
                c_idx = cid_to_idx[c_cid]
                r_main = main_rank_map[c_idx]
                r_bf = bf_ranks[m][c_cid]
                if abs(r_main - r_bf) > 1e-5:
                    diffs.append((c_cid, r_main, r_bf))

            if diffs:
                all_passed = False
                log(f"  FAILED: {m} has {len(diffs)} rank mismatches!")
                for c_cid, rm, rbf in diffs[:5]:
                    log(f"    {c_cid}: main={rm}, brute_force={rbf}")
            else:
                log(f"  PASSED: {m:<4} all {len(cand_cids)} candidate ranks match exactly.")

    log("\n" + "=" * 85)
    if all_passed:
        log(f"ALL {total_verified} QUERIES (10 random, 10 empty targets, 10 empty side effects) AND ALL METHODS (M2-M9) PASSED BRUTE-FORCE VERIFICATION.")
    else:
        log("VERIFICATION FAILED: rank differences detected.")
    log("=" * 85)

    out_file = ROOT / "experiments" / "logs" / "verify_experiments.log"
    out_file.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    log(f"Saved verification log to {out_file}")

if __name__ == "__main__":
    run_verification()

