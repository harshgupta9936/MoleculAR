"""
experiments/step4_weight_sensitivity.py

Step 4: Weight Sensitivity Grid Search on P2a (theta=0.4)
- 50/50 tune/test split of query drugs with seed 42.
- 66 triples (w_s, w_t, w_e) with step 0.1 summing to 1.0.
- Evaluates on tune set targeting MRR.
- Evaluates best triple vs baseline (0.4, 0.4, 0.2) on test set (MRR & enrichment with 95% bootstrap CIs).
- Saves experiments/results/weight_grid.csv.
- Saves experiments/figures/weight_heatmap.png.
"""

import csv
import sys
import time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_atc_l3, precompute_matrices,
    compute_average_ranks, rand_hit_at_k, bootstrap_ci, bootstrap_enrichment,
    DB_A_PATH, GOLD_PRIMARY_PATH, RESULTS_DIR
)

FIGURES_DIR = ROOT / "experiments" / "figures"

def generate_weight_triples():
    triples = []
    for i in range(11):
        for j in range(11 - i):
            k = 10 - i - j
            w_s = round(i / 10.0, 1)
            w_t = round(j / 10.0, 1)
            w_e = round(k / 10.0, 1)
            triples.append((w_s, w_t, w_e))
    return triples

def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== STEP 4: WEIGHT SENSITIVITY GRID SEARCH ===")
    t0 = time.time()

    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    matrices = precompute_matrices(drugs, pool_cids)

    cid_to_idx = matrices["cid_to_idx"]
    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    penalties = matrices["penalties"]
    N = len(pool_cids)

    # 1. Build P2a (0.4) queries
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}

    valid_queries = []
    for q_cid in cids_with_inds:
        q_idx = cid_to_idx[q_cid]
        q_inds = cid_ind_sets[q_cid]
        rel_P1 = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(cid_ind_sets[c])]
        rel_P2a = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4]
        if not rel_P2a:
            continue
        rel_indices = np.array([cid_to_idx[c] for c in rel_P2a], dtype=np.int32)
        cand_indices = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)
        valid_queries.append({
            "q_cid": q_cid,
            "q_idx": q_idx,
            "rel_indices": rel_indices,
            "cand_indices": cand_indices,
            "M": len(cand_indices),
            "R": len(rel_indices),
        })

    print(f"Total P2a (0.4) valid queries: {len(valid_queries)}")

    # 2. 50/50 tune/test split with seed 42 (split by drug)
    rng = np.random.RandomState(42)
    query_cids = [q["q_cid"] for q in valid_queries]
    perm = rng.permutation(len(query_cids))
    half = len(query_cids) // 2

    tune_indices = perm[:half]
    test_indices = perm[half:]

    tune_queries = [valid_queries[i] for i in tune_indices]
    test_queries = [valid_queries[i] for i in test_indices]

    print(f"Tune queries: {len(tune_queries)}, Test queries: {len(test_queries)}")

    # 3. Grid Search on Tune Set (66 triples)
    triples = generate_weight_triples()
    assert len(triples) == 66, f"Expected 66 triples, got {len(triples)}"

    grid_results = []
    best_triple = None
    best_mrr = -1.0

    print("Evaluating 66 weight triples on tune set...")
    for idx, (ws, wt, we) in enumerate(triples, 1):
        mrrs = []
        for q in tune_queries:
            q_idx = q["q_idx"]
            rel_idx = q["rel_indices"]
            cand_idx = q["cand_indices"]

            raw_scores = ws * S[q_idx] + wt * T[q_idx] + we * E[q_idx]
            final_scores = raw_scores / penalties

            rank_map = compute_average_ranks(final_scores, cand_idx)
            r_min = np.min([rank_map[i] for i in rel_idx])
            mrrs.append(1.0 / r_min)

        mean_mrr = float(np.mean(mrrs))
        grid_results.append({
            "w_s": ws,
            "w_t": wt,
            "w_e": we,
            "mrr_tune": mean_mrr,
        })

        if mean_mrr > best_mrr:
            best_mrr = mean_mrr
            best_triple = (ws, wt, we)

    print(f"Grid search complete. Best triple on tune set: w_s={best_triple[0]}, w_t={best_triple[1]}, w_e={best_triple[2]} with MRR = {best_mrr:.6f}")

    # Save weight_grid.csv
    grid_csv = RESULTS_DIR / "weight_grid.csv"
    with open(grid_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["w_s", "w_t", "w_e", "mrr_tune"])
        for r in grid_results:
            w.writerow([f"{r['w_s']:.1f}", f"{r['w_t']:.1f}", f"{r['w_e']:.1f}", f"{r['mrr_tune']:.6f}"])
    print(f"Saved weight grid to {grid_csv}")

    # 4. Generate Heatmap
    heatmap_matrix = np.full((11, 11), np.nan)
    for r in grid_results:
        # i = int(round(r["w_t"] * 10)) -> y axis
        # j = int(round(r["w_s"] * 10)) -> x axis
        i = int(round(r["w_t"] * 10))
        j = int(round(r["w_s"] * 10))
        heatmap_matrix[i, j] = r["mrr_tune"]

    plt.figure(figsize=(8, 7), dpi=300)
    cmap = plt.cm.viridis
    cmap.set_bad(color="lightgray")
    im = plt.imshow(heatmap_matrix, origin="lower", extent=[-0.05, 1.05, -0.05, 1.05], cmap=cmap)
    cbar = plt.colorbar(im)
    cbar.set_label("Mean MRR (Tune Set)", fontsize=11)

    plt.scatter([best_triple[0]], [best_triple[1]], color="red", marker="*", s=200, label=f"Best ({best_triple[0]}, {best_triple[1]}, {best_triple[2]}): {best_mrr:.4f}")
    plt.scatter([0.4], [0.4], color="cyan", marker="o", s=80, edgecolors="black", label=f"Default (0.4, 0.4, 0.2)")

    plt.title("Weight Sensitivity Grid Search on P2a ($\\theta=0.4$, Tune Set)", fontsize=12, fontweight="bold")
    plt.xlabel("Structural Weight ($w_s$)", fontsize=11)
    plt.ylabel("Target Weight ($w_t$)", fontsize=11)
    plt.xticks(np.arange(0.0, 1.1, 0.1))
    plt.yticks(np.arange(0.0, 1.1, 0.1))
    plt.grid(color="white", linestyle="--", linewidth=0.5, alpha=0.5)
    plt.legend(loc="upper right", framealpha=0.9)
    plt.tight_layout()

    heatmap_path = FIGURES_DIR / "weight_heatmap.png"
    plt.savefig(heatmap_path, dpi=300)
    plt.close()
    print(f"Saved heatmap to {heatmap_path}")

    # 5. Evaluate Best Triple vs Baseline (0.4, 0.4, 0.2) on Test Set
    print("\n--- Evaluating Best Triple vs Baseline on Test Set ---")
    triples_to_eval = [
        ("Best Tune Triple", best_triple),
        ("Default Baseline", (0.4, 0.4, 0.2)),
    ]

    test_comparison_rows = []
    for label, (ws, wt, we) in triples_to_eval:
        hits10 = []
        mrrs = []
        rand_hits10 = []
        for q in test_queries:
            q_idx = q["q_idx"]
            rel_idx = q["rel_indices"]
            cand_idx = q["cand_indices"]

            raw_scores = ws * S[q_idx] + wt * T[q_idx] + we * E[q_idx]
            final_scores = raw_scores / penalties

            rank_map = compute_average_ranks(final_scores, cand_idx)
            r_min = np.min([rank_map[i] for i in rel_idx])

            h10 = 1.0 if r_min <= 10.0 else 0.0
            mrr = 1.0 / r_min
            rnd = rand_hit_at_k(q["M"], q["R"], 10)

            hits10.append(h10)
            mrrs.append(mrr)
            rand_hits10.append(rnd)

        hits10 = np.array(hits10)
        mrrs = np.array(mrrs)
        rand_hits10 = np.array(rand_hits10)

        m_mrr, l_mrr, h_mrr = bootstrap_ci(mrrs, n_boot=1000, seed=42)
        m_enr, l_enr, h_enr = bootstrap_enrichment(hits10, rand_hits10, n_boot=1000, seed=42)

        print(f"{label} ({ws:.1f}, {wt:.1f}, {we:.1f}) on Test Set (N={len(test_queries)}):")
        print(f"  MRR:        {m_mrr:.4f} [{l_mrr:.4f}, {h_mrr:.4f}]")
        print(f"  Enrichment: {m_enr:.4f} [{l_enr:.4f}, {h_enr:.4f}]")

        test_comparison_rows.append({
            "model": label, "w_s": ws, "w_t": wt, "w_e": we, "n_test": len(test_queries),
            "mrr_mean": m_mrr, "mrr_low": l_mrr, "mrr_high": h_mrr,
            "enrich_mean": m_enr, "enrich_low": l_enr, "enrich_high": h_enr,
        })

    # Save test comparison to CSV
    test_csv = RESULTS_DIR / "weight_test_comparison.csv"
    with open(test_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["configuration", "w_s", "w_t", "w_e", "n_test", "mrr_mean", "mrr_low", "mrr_high", "enrich_mean", "enrich_low", "enrich_high"])
        for r in test_comparison_rows:
            w.writerow([
                r["model"], f"{r['w_s']:.1f}", f"{r['w_t']:.1f}", f"{r['w_e']:.1f}", r["n_test"],
                f"{r['mrr_mean']:.6f}", f"{r['mrr_low']:.6f}", f"{r['mrr_high']:.6f}",
                f"{r['enrich_mean']:.6f}", f"{r['enrich_low']:.6f}", f"{r['enrich_high']:.6f}"
            ])
    print(f"Saved test comparison to {test_csv}")
    print(f"Step 4 completed in {time.time() - t0:.2f}s.")

if __name__ == "__main__":
    main()

