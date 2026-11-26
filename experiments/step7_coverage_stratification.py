"""
experiments/step7_coverage_stratification.py

Executes Requirement N3:
Coverage-stratified evaluation for:
- P1
- P2a(0.4)
- P2b(0.4)
- P3-dn

Stratifies queries into four mutually exclusive groups:
  (a) has targets and has side effects
  (b) has targets only
  (c) has side effects only
  (d) has neither

Evaluates Methods M2 through M7 on:
  - Hit@10 (95% bootstrap CI)
  - MRR (95% bootstrap CI)
  - Enrichment (95% bootstrap CI)

Output saved to:
  - experiments/results/coverage_stratified_P1.csv
  - experiments/results/coverage_stratified_P2a_0.4.csv
  - experiments/results/coverage_stratified_P2b_0.4.csv
  - experiments/results/coverage_stratified_P3_dn.csv
  - experiments/results/coverage_stratified_all.csv
"""

import csv
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_atc_l3, precompute_matrices,
    compute_average_ranks, bootstrap_ci, bootstrap_enrichment,
    rand_hit_at_k, load_duplicates,
    DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

def evaluate_coverage_stratification():
    print("="*80)
    print("STEP 7: COVERAGE-STRATIFIED EVALUATION (N3)")
    print("="*80)
    t0 = time.time()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load data
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    cid_atc_l3 = load_atc_l3(set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids, DB_A_PATH)
    matrices = precompute_matrices(drugs, pool_cids)

    cid_to_idx = matrices["cid_to_idx"]
    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    penalties = matrices["penalties"]
    N = len(pool_cids)

    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}

    pool_with_atc = {cid for cid in pool_cids if len(cid_atc_l3[cid]) > 0}
    atc_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc], dtype=np.int32)

    protocols = ["P1", "P2a_0.4", "P2b_0.4", "P3_dn"]
    methods = ["M2", "M3", "M4", "M5", "M6", "M7"]

    all_stratified_rows = []

    for proto in protocols:
        print(f"\n--- Processing Protocol: {proto} ---")
        queries = []

        for q_cid in cids_with_inds:
            q_idx = cid_to_idx[q_cid]
            q_inds = cid_ind_sets[q_cid]

            # Base P1 relevant
            rel_P1 = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(cid_ind_sets[c])]
            if not rel_P1:
                continue

            if proto == "P1":
                rel_cids = rel_P1
                cand_indices = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)

            elif proto == "P2a_0.4":
                rel_cids = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4]
                if not rel_cids:
                    continue
                cand_indices = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)

            elif proto == "P2b_0.4":
                rel_cids = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4]
                if not rel_cids:
                    continue
                cand_indices = np.array([i for i in range(N) if i != q_idx and S[q_idx, i] <= 0.4], dtype=np.int32)

            elif proto == "P3_dn":
                if q_cid not in pool_with_atc:
                    continue
                q_atc = cid_atc_l3[q_cid]
                dups_q = duplicates_of.get(q_cid, set())
                rel_cids = [c for c in rel_P1 if c in pool_with_atc and q_atc.isdisjoint(cid_atc_l3[c]) and c not in dups_q]
                if not rel_cids:
                    continue
                cand_indices = np.array([i for i in atc_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)

            else:
                raise ValueError(f"Unknown protocol {proto}")

            has_t = len(drugs[q_cid]["targets"]) > 0
            has_se = len(drugs[q_cid]["side_effects"]) > 0

            if has_t and has_se:
                group = "has_targets_and_se"
                group_label = "(a) Targets and Side Effects"
            elif has_t and not has_se:
                group = "targets_only"
                group_label = "(b) Targets Only"
            elif not has_t and has_se:
                group = "se_only"
                group_label = "(c) Side Effects Only"
            else:
                group = "neither"
                group_label = "(d) Neither"

            queries.append({
                "q_cid": q_cid,
                "q_idx": q_idx,
                "cand_indices": cand_indices,
                "rel_indices": np.array([cid_to_idx[c] for c in rel_cids], dtype=np.int32),
                "M": len(cand_indices),
                "R": len(rel_cids),
                "group": group,
                "group_label": group_label,
            })

        print(f"Total queries for {proto}: {len(queries)}")

        # Compute query metrics for each method
        query_metrics = {m: {"hit10": [], "mrr": [], "rand_hit10": []} for m in methods}
        query_groups = [q["group"] for q in queries]

        for q in queries:
            q_idx = q["q_idx"]
            cand_indices = q["cand_indices"]
            rel_indices = q["rel_indices"]
            M = q["M"]
            R = q["R"]
            rand_h10 = rand_hit_at_k(M, R, 10)

            s_row = S[q_idx]
            t_row = T[q_idx]
            e_row = E[q_idx]

            scores_m = {
                "M2": s_row,
                "M3": t_row,
                "M4": e_row,
                "M5": 0.5 * s_row + 0.5 * t_row,
                "M6": 0.4 * s_row + 0.4 * t_row + 0.2 * e_row,
                "M7": (0.4 * s_row + 0.4 * t_row + 0.2 * e_row) / penalties,
            }

            for m in methods:
                rank_map = compute_average_ranks(scores_m[m], cand_indices)
                rel_ranks = [rank_map[idx] for idx in rel_indices]
                r_min = min(rel_ranks)
                hit10 = 1.0 if r_min <= 10.0 else 0.0
                mrr = 1.0 / r_min

                query_metrics[m]["hit10"].append(hit10)
                query_metrics[m]["mrr"].append(mrr)
                query_metrics[m]["rand_hit10"].append(rand_h10)

        # Convert to numpy arrays
        for m in methods:
            query_metrics[m]["hit10"] = np.array(query_metrics[m]["hit10"], dtype=np.float64)
            query_metrics[m]["mrr"] = np.array(query_metrics[m]["mrr"], dtype=np.float64)
            query_metrics[m]["rand_hit10"] = np.array(query_metrics[m]["rand_hit10"], dtype=np.float64)

        groups_ordered = [
            ("has_targets_and_se", "(a) Targets and Side Effects"),
            ("targets_only", "(b) Targets Only"),
            ("se_only", "(c) Side Effects Only"),
            ("neither", "(d) Neither"),
        ]

        proto_rows = []

        for grp_key, grp_name in groups_ordered:
            sub_indices = [i for i, g in enumerate(query_groups) if g == grp_key]
            n_sub = len(sub_indices)
            print(f"  Group {grp_name}: n = {n_sub}")
            if n_sub == 0:
                continue

            sub_idx_arr = np.array(sub_indices, dtype=np.int32)

            for m in methods:
                h10_sub = query_metrics[m]["hit10"][sub_idx_arr]
                mrr_sub = query_metrics[m]["mrr"][sub_idx_arr]
                rand_h10_sub = query_metrics[m]["rand_hit10"][sub_idx_arr]

                h10_mean, h10_low, h10_high = bootstrap_ci(h10_sub, n_boot=1000, seed=42)
                mrr_mean, mrr_low, mrr_high = bootstrap_ci(mrr_sub, n_boot=1000, seed=42)
                enr_mean, enr_low, enr_high = bootstrap_enrichment(h10_sub, rand_h10_sub, n_boot=1000, seed=42)

                row_base = {
                    "protocol": proto,
                    "group_code": grp_key,
                    "group_name": grp_name,
                    "n_queries": n_sub,
                    "method": m,
                }

                r_h10 = dict(row_base, metric="hit@10", mean=h10_mean, ci_low=h10_low, ci_high=h10_high)
                r_mrr = dict(row_base, metric="mrr", mean=mrr_mean, ci_low=mrr_low, ci_high=mrr_high)
                r_enr = dict(row_base, metric="enrichment", mean=enr_mean, ci_low=enr_low, ci_high=enr_high)

                proto_rows.extend([r_h10, r_mrr, r_enr])
                all_stratified_rows.extend([r_h10, r_mrr, r_enr])

        # Save per-protocol CSV
        proto_csv = RESULTS_DIR / f"coverage_stratified_{proto}.csv"
        with open(proto_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["protocol", "group_code", "group_name", "n_queries", "method", "metric", "mean", "ci_low", "ci_high"])
            w.writeheader()
            w.writerows(proto_rows)
        print(f"Saved {proto_csv} ({len(proto_rows)} rows)")

    # Save consolidated CSV
    all_csv = RESULTS_DIR / "coverage_stratified_all.csv"
    with open(all_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["protocol", "group_code", "group_name", "n_queries", "method", "metric", "mean", "ci_low", "ci_high"])
        w.writeheader()
        w.writerows(all_stratified_rows)
    print(f"\nSaved consolidated {all_csv} ({len(all_stratified_rows)} rows) in {time.time() - t0:.2f}s")

if __name__ == "__main__":
    evaluate_coverage_stratification()
