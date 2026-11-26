"""
experiments/generate_findings_tables.py

Generates experiments/results/FINDINGS_TABLES.md strictly from on-disk CSV files:
- experiments/results/paired_differences.csv
- experiments/results/weight_stability_p2b_dn.csv
- experiments/results/weight_stability_p1_dn.csv
- experiments/results/hub_penalty_ablation.csv
- experiments/results/failure_analysis_p2b_dn.csv
- experiments/results/failure_analysis_p3_dn.csv

Under every table, prints the source CSV path and SHA-256 hash.
"""

import csv
import hashlib
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "experiments" / "results"
OUT_FILE = RESULTS_DIR / "FINDINGS_TABLES.md"

def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    lines = []
    def p(l=""):
        lines.append(l)

    p("# MoleculAR Second-Round Findings & Analyses (Duplicate-Neutral)")
    p()
    p("> **Generated strictly from on-disk CSV files by:** `experiments/generate_findings_tables.py`  ")
    p("> **Compliance:** Fulfills Requirements N1 and N2 (no hand-typed numbers; source CSV & SHA-256 printed under every table).")
    p()
    p("---")
    p()

    # =========================================================================
    # 1. Paired Comparisons
    # =========================================================================
    p("## 1. Paired Bootstrap Comparisons (M6 vs M5, M7 vs M6, M7 vs M5, M9 vs M7, M7 vs M2, M7 vs M8, M6 vs M2)")
    p()
    p("Paired bootstrap differences evaluated on Configuration C1 across duplicate-neutral protocols (`P1-dn`, `P2a-dn (0.4)`, `P2b-dn (0.4)`, `P3-dn`) with 1,000 paired resamples (seed 42).")
    p()

    paired_csv = RESULTS_DIR / "paired_differences.csv"
    paired_sha = sha256_file(paired_csv)

    # Load paired diffs
    paired_data = []
    with open(paired_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            paired_data.append(r)

    # Group by protocol
    protocols = ["P1-dn", "P2a-dn (0.4)", "P2b-dn (0.4)", "P3-dn"]
    for proto in protocols:
        p(f"### Table 1.{proto}: Paired Differences on {proto}")
        p()
        p("| Comparison | Metric | Mean Difference | 95% Bootstrap CI | Significant (CI excludes 0)? |")
        p("| :--- | :---: | :---: | :---: | :---: |")
        proto_rows = [r for r in paired_data if r["protocol"] == proto]
        for r in proto_rows:
            mean_d = float(r["mean_diff"])
            ci_l = float(r["ci_low"])
            ci_h = float(r["ci_high"])
            sig = "Yes" if (ci_l > 0 or ci_h < 0) else "No"
            metric_label = "MRR" if r["metric"] == "mrr" else "Hit@10"
            p(f"| {r['comparison']} | {metric_label} | {mean_d:+.4f} | [{ci_l:+.4f}, {ci_h:+.4f}] | {sig} |")
        p()
        p(f"> **Source File:** [`{paired_csv.relative_to(ROOT)}`](file:///{paired_csv.resolve().as_posix()})  ")
        p(f"> **SHA-256:** `{paired_sha}`")
        p()

    p("---")
    p()

    # =========================================================================
    # 2. Weight Stability Analysis
    # =========================================================================
    p("## 2. Weight Stability Analysis (20 Random 50/50 Splits, Seeds 0–19)")
    p()
    p("Evaluates whether optimal modality weights $(w_S, w_T, w_E)$ selected on a 50% tuning split over the 66-simplex grid generalize to an independent 50% test split, compared against production default weights $(0.4, 0.4, 0.2)$.")
    p()

    for p_name, csv_name in [("P2b-dn (0.4)", "weight_stability_p2b_dn.csv"), ("P1-dn", "weight_stability_p1_dn.csv")]:
        ws_csv = RESULTS_DIR / csv_name
        ws_sha = sha256_file(ws_csv)

        ws_rows = []
        with open(ws_csv, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                ws_rows.append(r)

        p(f"### Table 2.{p_name}: Weight Stability on {p_name}")
        p()
        p("| Split Seed | Best Tune Weights $(w_S, w_T, w_E)$ | Best Tune MRR | Best Test MRR | Default Test MRR | Test MRR Diff (Best - Def) | Default Tune MRR | Default Tune Rank (/66) |")
        p("| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")

        diffs_mrr = []
        best_test_mrrs = []
        def_test_mrrs = []
        for r in ws_rows:
            seed = int(r["split_seed"])
            w_str = f"({float(r['best_w_s']):.1f}, {float(r['best_w_t']):.1f}, {float(r['best_w_e']):.1f})"
            b_tune = float(r["best_tune_mrr"])
            b_test = float(r["best_test_mrr"])
            d_test = float(r["default_test_mrr"])
            d_tune = float(r["default_tune_mrr"])
            d_diff = float(r["test_diff_mrr"])
            d_rank = int(float(r["default_tune_rank"]))
            diffs_mrr.append(d_diff)
            best_test_mrrs.append(b_test)
            def_test_mrrs.append(d_test)
            p(f"| {seed} | {w_str} | {b_tune:.4f} | {b_test:.4f} | {d_test:.4f} | {d_diff:+.4f} | {d_tune:.4f} | {d_rank} |")

        mean_diff = np.mean(diffs_mrr)
        std_diff = np.std(diffs_mrr)
        mean_best = np.mean(best_test_mrrs)
        mean_def = np.mean(def_test_mrrs)

        p(f"| **Mean** | — | — | **{mean_best:.4f}** | **{mean_def:.4f}** | **{mean_diff:+.4f} (std: {std_diff:.4f})** | — | — |")
        p()
        p(f"> **Source File:** [`{ws_csv.relative_to(ROOT)}`](file:///{ws_csv.resolve().as_posix()})  ")
        p(f"> **SHA-256:** `{ws_sha}`")
        p()

    p("---")
    p()

    # =========================================================================
    # 3. Hub Penalty Exponent Ablation & Recurrence Statistics
    # =========================================================================
    p("## 3. Hub Penalty Exponent Ablation and Hubness Recurrence Statistics")
    p()
    p("Evaluates the impact of hub penalty exponent $k \\in \\{0, 1, 2, 3, 4\\}$ on ranking performance and hub drug dominance.")
    p()

    hub_csv = RESULTS_DIR / "hub_penalty_ablation.csv"
    hub_sha = sha256_file(hub_csv)

    hub_rows = []
    with open(hub_csv, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            hub_rows.append(r)

    p("### Table 3.1: Hub Penalty Exponent Ablation across Duplicate-Neutral Protocols")
    p()
    p("| Protocol | $k$ | N | MRR (95% CI) | Hit@10 (95% CI) | Enrichment (95% CI) | Mean Top-10 Hubness | Max Candidate Appearances | Max Candidate Share | Top 1% Candidate Share |")
    p("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")

    for r in hub_rows:
        proto = r["protocol"]
        k = int(r["k"])
        n_q = int(r["n_queries"])
        mrr_str = f"{float(r['mrr']):.4f} [{float(r['mrr_ci_low']):.4f}, {float(r['mrr_ci_high']):.4f}]"
        h10_str = f"{float(r['hit10']):.4f} [{float(r['hit10_ci_low']):.4f}, {float(r['hit10_ci_high']):.4f}]"
        enr_str = f"{float(r['enrichment']):.2f} [{float(r['enrichment_ci_low']):.2f}, {float(r['enrichment_ci_high']):.2f}]"
        hub_mean = f"{float(r['mean_top10_hubness']):.2f}"
        max_cnt = int(r["max_cand_count"])
        max_sh = f"{float(r['max_share']) * 100:.2f}%"
        top1_sh = f"{float(r['top1pct_share']) * 100:.2f}%"
        p(f"| {proto} | {k} | {n_q} | {mrr_str} | {h10_str} | {enr_str} | {hub_mean} | {max_cnt} | {max_sh} | {top1_sh} |")

    p()
    p(f"> **Source File:** [`{hub_csv.relative_to(ROOT)}`](file:///{hub_csv.resolve().as_posix()})  ")
    p(f"> **SHA-256:** `{hub_sha}`")
    p()
    p("---")
    p()

    # =========================================================================
    # 4. Failure Analysis
    # =========================================================================
    p("## 4. Failure Analysis (Top 20 Best and Worst Queries for Method M7)")
    p()

    for p_name, csv_name in [("P2b-dn (0.4)", "failure_analysis_p2b_dn.csv"), ("P3-dn", "failure_analysis_p3_dn.csv")]:
        fa_csv = RESULTS_DIR / csv_name
        fa_sha = sha256_file(fa_csv)

        fa_rows = []
        with open(fa_csv, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                fa_rows.append(r)

        best_rows = [r for r in fa_rows if r["category"] == "best_20"]
        worst_rows = [r for r in fa_rows if r["category"] == "worst_20"]

        p(f"### Table 4.{p_name}.1: Top 20 Best Performing Queries for M7 on {p_name}")
        p()
        p("| Query Drug | N (Rel) | First Rel Drug | Rank ($r_{\\min}$) | $S$ | $T$ | $E$ | Targets (Q / Cand)? | SEs (Q / Cand)? |")
        p("| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for r in best_rows:
            q_str = f"{r['query_cid']} ({r['query_name']})"
            rel_str = f"{r['first_rel_cid']} ({r['first_rel_name']})"
            rank_str = f"{float(r['first_rel_rank']):.1f}"
            s_val = f"{float(r['S']):.3f}"
            t_val = f"{float(r['T']):.3f}"
            e_val = f"{float(r['E']):.3f}"
            t_flag = f"{'Y' if r['query_has_targets'] == 'True' else 'N'} / {'Y' if r['cand_has_targets'] == 'True' else 'N'}"
            se_flag = f"{'Y' if r['query_has_se'] == 'True' else 'N'} / {'Y' if r['cand_has_se'] == 'True' else 'N'}"
            p(f"| {q_str} | {r['n_relevant']} | {rel_str} | {rank_str} | {s_val} | {t_val} | {e_val} | {t_flag} | {se_flag} |")

        p()
        p(f"> **Source File:** [`{fa_csv.relative_to(ROOT)}`](file:///{fa_csv.resolve().as_posix()})  ")
        p(f"> **SHA-256:** `{fa_sha}`")
        p()

        p(f"### Table 4.{p_name}.2: Top 20 Worst Performing Queries for M7 on {p_name}")
        p()
        p("| Query Drug | N (Rel) | First Rel Drug | Rank ($r_{\\min}$) | $S$ | $T$ | $E$ | Targets (Q / Cand)? | SEs (Q / Cand)? |")
        p("| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for r in worst_rows:
            q_str = f"{r['query_cid']} ({r['query_name']})"
            rel_str = f"{r['first_rel_cid']} ({r['first_rel_name']})"
            rank_str = f"{float(r['first_rel_rank']):.1f}"
            s_val = f"{float(r['S']):.3f}"
            t_val = f"{float(r['T']):.3f}"
            e_val = f"{float(r['E']):.3f}"
            t_flag = f"{'Y' if r['query_has_targets'] == 'True' else 'N'} / {'Y' if r['cand_has_targets'] == 'True' else 'N'}"
            se_flag = f"{'Y' if r['query_has_se'] == 'True' else 'N'} / {'Y' if r['cand_has_se'] == 'True' else 'N'}"
            p(f"| {q_str} | {r['n_relevant']} | {rel_str} | {rank_str} | {s_val} | {t_val} | {e_val} | {t_flag} | {se_flag} |")

        p()
        p(f"> **Source File:** [`{fa_csv.relative_to(ROOT)}`](file:///{fa_csv.resolve().as_posix()})  ")
        p(f"> **SHA-256:** `{fa_sha}`")
        p()

    OUT_FILE.write_text("\n".join(lines), encoding="utf-8")
    print(f"Generated {OUT_FILE} ({len(lines)} lines)")

if __name__ == "__main__":
    main()
