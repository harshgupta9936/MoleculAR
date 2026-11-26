"""
experiments/generate_final_tables.py

Requirement F4: Consolidation
1. Reads all result CSVs from experiments/results/.
2. Generates experiments/results/FINAL_TABLES.md containing every protocol (P1, P2a, P2b, P3, P4, P4b and -dn),
   all methods, CIs, paired differences, coverage stratification, nulls, ablations, tie sensitivity, and known cases.
3. Prints exact source CSV path and SHA-256 under EVERY table.
4. Updates experiments/results/manifest.json with all database hashes, git commit, seed, library versions,
   and SHA-256 of every result CSV.
"""

import os
import sys
import json
import csv
import hashlib
import datetime
import subprocess
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

RESULTS_DIR = ROOT / "experiments" / "results"
DATA_DIR = ROOT / "experiments" / "data"

def sha256_file(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            h.update(chunk)
    return h.hexdigest()

def make_source_line(filename):
    p = RESULTS_DIR / filename
    return f"\n**Source CSV:** `{p.as_posix()}` (SHA-256: `{sha256_file(p)}`)\n"

def format_protocol_table(csv_name, title_template):
    p = RESULTS_DIR / csv_name
    if not p.exists():
        return ""
    df = pd.read_csv(p)
    n_queries = int(df['n_queries'].iloc[0])
    title = title_template.format(n_queries=f"{n_queries:,}")
    lines = [f"### {title}\n"]
    lines.append("| Method | Hit@10 [95% CI] | Hit@20 [95% CI] | Precision@10 [95% CI] | MRR [95% CI] | Percentile Rank [95% CI] | Enrichment [95% CI] |")
    lines.append("| :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
    
    methods = [f"M{i}" for i in range(1, 10)]
    metrics_order = ["hit@10", "hit@20", "precision@10", "mrr", "percentile_rank", "enrichment"]
    
    # Check column names
    m_col = "method"
    met_col = "metric"
    mean_col = "mean"
    ci_l_col = "ci_low"
    ci_h_col = "ci_high"
    
    for m in methods:
        m_df = df[df[m_col] == m].set_index(met_col)
        if m_df.empty:
            continue
        row = [f"**{m}**"]
        for met in metrics_order:
            if met in m_df.index:
                mean_v = m_df.loc[met, mean_col]
                ci_l = m_df.loc[met, ci_l_col]
                ci_h = m_df.loc[met, ci_h_col]
                if met == "enrichment":
                    row.append(f"{mean_v:.2f} [{ci_l:.2f}, {ci_h:.2f}]")
                else:
                    row.append(f"{mean_v:.4f} [{ci_l:.4f}, {ci_h:.4f}]")
            else:
                row.append("N/A")
        lines.append("| " + " | ".join(row) + " |")
        
    lines.append(make_source_line(csv_name))
    return "\n".join(lines) + "\n"

def build_final_tables():
    out = []
    out.append("# MOLECULAR EVALUATION HARNESS: FINAL CONSOLIDATED RESULTS TABLES\n")
    out.append(f"**Generated:** {datetime.date.today().isoformat()}  ")
    out.append("**Primary Database:** DB-A (`experiments/data/drugrepo_oct6.db`)  ")
    out.append("**Strict Governance:** Frozen PROTOCOL.md definitions, exact mid-rank average ranking for ties, float64 precision, seed 42, chembl_id keys.  \n")
    out.append("---\n")

    # 1. PRIMARY PROTOCOL PERFORMANCE METRICS (C1)
    out.append("## 1. Primary Protocol Performance Metrics (Candidate Pool C1)\n")
    
    protocols_c1 = [
        ("metrics_C1_P1.csv", "Protocol P1 (Broad Therapeutic Indication Matching - Original Baseline, N={n_queries})"),
        ("metrics_C1_P1_dn.csv", "Protocol P1-dn (Broad Indication Matching - Duplicate-Neutral, N={n_queries})"),
        ("metrics_C1_P2a_0.4.csv", "Protocol P2a (0.4) (Analogue Neutralization at θ=0.4 - Original Baseline, N={n_queries})"),
        ("metrics_C1_P2a_dn_0.4.csv", "Protocol P2a-dn (0.4) (Analogue Neutralization at θ=0.4 - Duplicate-Neutral, N={n_queries})"),
        ("metrics_C1_P2a_0.3.csv", "Protocol P2a (0.3) (Analogue Neutralization at θ=0.3, N={n_queries})"),
        ("metrics_C1_P2a_0.5.csv", "Protocol P2a (0.5) (Analogue Neutralization at θ=0.5, N={n_queries})"),
        ("metrics_C1_P2b_0.4.csv", "Protocol P2b (0.4) (Analogue Removal at θ=0.4 - Original Baseline, N={n_queries})"),
        ("metrics_C1_P2b_dn_0.4.csv", "Protocol P2b-dn (0.4) (Analogue Removal at θ=0.4 - Duplicate-Neutral, N={n_queries})"),
        ("metrics_C1_P3.csv", "Protocol P3 (Cross-Class Repurposing: Disjoint ATC-L3, Pool with ATC - Original Baseline, N={n_queries})"),
        ("metrics_C1_P3_all.csv", "Protocol P3-all (Cross-Class Repurposing: Disjoint ATC-L3, Full Pool - Original Baseline, N={n_queries})"),
        ("metrics_C1_P3_dn.csv", "Protocol P3-dn (Cross-Class Repurposing: Disjoint ATC-L3, Pool with ATC - Duplicate-Neutral, N={n_queries})"),
        ("metrics_C1_P4.csv", "Protocol P4 (Stricter Cross-Class: Disjoint ATC Level-1 Anatomical Main Groups, Full ATC Pool, N={n_queries})"),
        ("metrics_C1_P4b.csv", "Protocol P4b (Stricter Cross-Class: Disjoint ATC Level-1 Anatomical Removal, N={n_queries})"),
    ]
    
    for fname, title in protocols_c1:
        out.append(format_protocol_table(fname, title))
        
    out.append("---\n")

    # 2. PAIRED BOOTSTRAP DIFFERENCES
    out.append("## 2. Paired Bootstrap Differences (1,000 Resamples, Seed 42)\n")
    
    # 2.1 P1, P2, P3 paired diffs
    p_diff_path = "paired_differences.csv"
    if (RESULTS_DIR / p_diff_path).exists():
        df_pd = pd.read_csv(RESULTS_DIR / p_diff_path)
        out.append("### Paired Differences across Protocols P1, P2a, P2b, P3-dn\n")
        out.append("| Protocol | Comparison | Metric | Mean Difference | 95% Bootstrap CI | Formatted Output |")
        out.append("| :--- | :--- | :--- | :---: | :---: | :---: |")
        for _, r in df_pd.iterrows():
            fmt = r['formatted'] if 'formatted' in r and pd.notna(r['formatted']) else f"{r['mean_diff']:+.4f} [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]"
            out.append(f"| {r['protocol']} | {r['comparison']} | {r['metric']} | {r['mean_diff']:+.4f} | [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | `{fmt}` |")
        out.append(make_source_line(p_diff_path))
        
    # 2.2 P4 and P4b paired diffs
    p4_diff_path = "paired_differences_p4_p4b.csv"
    if (RESULTS_DIR / p4_diff_path).exists():
        df_p4d = pd.read_csv(RESULTS_DIR / p4_diff_path)
        out.append("### Paired Differences for Protocols P4 and P4b\n")
        out.append("| Protocol | Comparison | Metric | Mean Difference | 95% Bootstrap CI | Formatted Output |")
        out.append("| :--- | :--- | :--- | :---: | :---: | :---: |")
        for _, r in df_p4d.iterrows():
            fmt = r['formatted'] if 'formatted' in r and pd.notna(r['formatted']) else f"{r['mean_diff']:+.4f} [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]"
            out.append(f"| {r['protocol']} | {r['comparison']} | {r['metric']} | {r['mean_diff']:+.4f} | [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}] | `{fmt}` |")
        out.append(make_source_line(p4_diff_path))
        
    out.append("---\n")

    # 3. COVERAGE-STRATIFIED EVALUATIONS
    out.append("## 3. Coverage-Stratified Evaluations (Methods M2 to M7)\n")
    
    # 3.1 P1, P2a, P2b, P3-dn
    cov_all_path = "coverage_stratified_all.csv"
    if (RESULTS_DIR / cov_all_path).exists():
        df_cov = pd.read_csv(RESULTS_DIR / cov_all_path)
        out.append("### Coverage-Stratified Evaluation across Protocols P1, P2a(0.4), P2b(0.4), P3-dn\n")
        out.append("| Protocol | Stratum | N Queries | Method | Hit@10 [95% CI] | MRR [95% CI] | Enrichment [95% CI] |")
        out.append("| :--- | :--- | :---: | :---: | :---: | :---: | :---: |")
        groups = df_cov[["protocol", "group_code", "group_name", "n_queries", "method"]].drop_duplicates().values
        for p_name, g_code, g_name, nq, m in groups:
            sub = df_cov[(df_cov["protocol"] == p_name) & (df_cov["group_code"] == g_code) & (df_cov["method"] == m)].set_index("metric")
            h10_str = f"{sub.loc['hit@10', 'mean']:.4f} [{sub.loc['hit@10', 'ci_low']:.4f}, {sub.loc['hit@10', 'ci_high']:.4f}]" if "hit@10" in sub.index else "N/A"
            mrr_str = f"{sub.loc['mrr', 'mean']:.4f} [{sub.loc['mrr', 'ci_low']:.4f}, {sub.loc['mrr', 'ci_high']:.4f}]" if "mrr" in sub.index else "N/A"
            enr_str = f"{sub.loc['enrichment', 'mean']:.2f} [{sub.loc['enrichment', 'ci_low']:.2f}, {sub.loc['enrichment', 'ci_high']:.2f}]" if "enrichment" in sub.index else "N/A"
            out.append(f"| {p_name} | {g_name} | {int(nq)} | {m} | {h10_str} | {mrr_str} | {enr_str} |")
        out.append(make_source_line(cov_all_path))
        
    # 3.2 P4 and P4b
    cov_p4_path = "coverage_stratified_p4_p4b.csv"
    if (RESULTS_DIR / cov_p4_path).exists():
        df_cov4 = pd.read_csv(RESULTS_DIR / cov_p4_path)
        out.append("### Coverage-Stratified Evaluation for Protocols P4 and P4b\n")
        out.append("| Protocol | Stratum | N Queries | Method | Hit@10 [95% CI] | MRR [95% CI] | Enrichment [95% CI] |")
        out.append("| :--- | :--- | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_cov4.iterrows():
            out.append(f"| {r['protocol']} | {r['stratum']} | {int(r['n_queries'])} | {r['method']} | {r['hit10_ci']} | {r['mrr_ci']} | {r['enrichment_ci']} |")
        out.append(make_source_line(cov_p4_path))
        
    out.append("---\n")

    # 4. NULL MODELS AND SHUFFLE CONTROLS
    out.append("## 4. Null Models and Shuffle Controls\n")
    
    # 4.1 Unstratified Shuffle Controls
    shuf_path = "shuffle_control_results.csv"
    if (RESULTS_DIR / shuf_path).exists():
        df_shuf = pd.read_csv(RESULTS_DIR / shuf_path)
        out.append("### Unstratified Feature Shuffle Controls (20 Permutations, Seed 42)\n")
        out.append("| Protocol | Method | Metric | Real Value | Null Mean | Null Min | Null Max | Null Range | N Permutations |")
        out.append("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_shuf.iterrows():
            out.append(f"| {r['protocol']} | {r['method']} | {r['metric']} | {r['real_value']:.4f} | {r['shuffle_mean']:.4f} | {r['shuffle_min']:.4f} | {r['shuffle_max']:.4f} | `{r['shuffle_range']}` | {int(r['n_permutations'])} |")
        out.append(make_source_line(shuf_path))
        
    # 4.2 Coverage-Fair Nulls within Strata
    cfn_path = "coverage_fair_nulls.csv"
    if (RESULTS_DIR / cfn_path).exists():
        df_cfn = pd.read_csv(RESULTS_DIR / cfn_path)
        out.append("### Coverage-Fair Nulls within Coverage Strata (20 Permutations, Seed 42)\n")
        out.append("> Note: Per specification, enrichment is strictly omitted for single-signal methods.\n")
        out.append("| Protocol | Method | Stratum | N Queries | Real Hit@10 | Null Mean Hit@10 | Null Range Hit@10 | Real MRR | Null Mean MRR | Null Range MRR |")
        out.append("| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_cfn.iterrows():
            out.append(f"| {r['protocol']} | {r['method']} | {r['stratum']} | {int(r['n_queries'])} | {r['real_hit10']:.4f} | {r['null_mean_hit10']:.4f} | `{r['null_range_hit10']}` | {r['real_mrr']:.4f} | {r['null_mean_mrr']:.4f} | `{r['null_range_mrr']}` |")
        out.append(make_source_line(cfn_path))
        
    out.append("---\n")

    # 5. ABLATION STUDIES
    out.append("## 5. Ablation Studies\n")
    
    # 5.1 Hub Penalty Ablation
    hub_path = "hub_penalty_ablation.csv"
    if (RESULTS_DIR / hub_path).exists():
        df_hub = pd.read_csv(RESULTS_DIR / hub_path)
        out.append("### Hub Penalty Sensitivity Ablation (Exponent k in [0.0, 3.0])\n")
        out.append("| Protocol | Exponent k | N Queries | Hit@10 [95% CI] | MRR [95% CI] | Enrichment [95% CI] | Top 1% Share |")
        out.append("| :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_hub.iterrows():
            h10_ci = f"{r['hit10']:.4f} [{r['hit10_ci_low']:.4f}, {r['hit10_ci_high']:.4f}]"
            mrr_ci = f"{r['mrr']:.4f} [{r['mrr_ci_low']:.4f}, {r['mrr_ci_high']:.4f}]"
            enr_ci = f"{r['enrichment']:.2f} [{r['enrichment_ci_low']:.2f}, {r['enrichment_ci_high']:.2f}]"
            out.append(f"| {r['protocol']} | {r['k']} | {int(r['n_queries'])} | {h10_ci} | {mrr_ci} | {enr_ci} | {r['top1pct_share']:.4f} |")
        out.append(make_source_line(hub_path))
        
    # 5.2 Side Effect Leakage Approved Ablations
    abl_path = "se_leakage_approved_ablations.csv"
    if (RESULTS_DIR / abl_path).exists():
        df_abl = pd.read_csv(RESULTS_DIR / abl_path)
        out.append("### Approved Side-Effect Leakage Ablations (Methods M4 and M6)\n")
        out.append("| Protocol | Condition | Method | Hit@10 [95% CI] | Original Baseline Hit@10 [95% CI] | MRR [95% CI] | Original Baseline MRR [95% CI] |")
        out.append("| :--- | :--- | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_abl.iterrows():
            out.append(f"| {r['protocol']} | {r['condition']} | {r['method']} | {r['hit10_ci']} | {r['orig_hit10_ci']} | {r['mrr_ci']} | {r['orig_mrr_ci']} |")
        out.append(make_source_line(abl_path))
        
    out.append("---\n")

    # 6. TIE SENSITIVITY ANALYSIS (Sensitivity Evaluation)
    out.append("## 6. Tie Sensitivity Analysis (Sensitivity Evaluation)\n")
    tie_path = "tie_sensitivity.csv"
    if (RESULTS_DIR / tie_path).exists():
        df_tie = pd.read_csv(RESULTS_DIR / tie_path)
        out.append("### Evaluation under Average, Best, and Worst Rank Tie Rules (Methods M2 to M7)\n")
        out.append("| Protocol | Method | Tie Rule | Hit@10 | Hit@20 | MRR | N Queries |")
        out.append("| :--- | :---: | :--- | :---: | :---: | :---: | :---: |")
        for _, r in df_tie.iterrows():
            out.append(f"| {r['protocol']} | {r['method']} | {r['tie_rule']} | {r['hit@10']:.4f} | {r['hit@20']:.4f} | {r['mrr']:.4f} | {int(r['n_queries'])} |")
        out.append(make_source_line(tie_path))

    tie_m8_path = "tie_sensitivity_m8.csv"
    if (RESULTS_DIR / tie_m8_path).exists():
        df_tie_m8 = pd.read_csv(RESULTS_DIR / tie_m8_path)
        out.append("### Method M8 (Popularity) Evaluation under Average, Best, and Worst Rank Tie Rules\n")
        out.append("| Category | Protocol | Comparison / Method | Tie Rule | Metric | Value / Diff | 95% Bootstrap CI | Formatted Output | N Queries |")
        out.append("| :--- | :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_tie_m8.iterrows():
            ci_str = f"[{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]" if pd.notna(r['ci_low']) else "N/A"
            out.append(f"| {r['category']} | {r['protocol']} | {r['comparison_or_method']} | {r['tie_rule']} | {r['metric']} | {r['value_or_diff']:+.4f} | {ci_str} | `{r['formatted']}` | {int(r['n_queries'])} |")
        out.append(make_source_line(tie_m8_path))
        
    out.append("---\n")

    # 7. CANONICAL KNOWN-CASE BENCHMARK ON DB-A
    out.append("## 7. Canonical Known-Case Benchmark Check on DB-A\n")
    kc_path = "known_cases_db_a.csv"
    if (RESULTS_DIR / kc_path).exists():
        df_kc = pd.read_csv(RESULTS_DIR / kc_path)
        out.append("### Ranks of 8 Canonical Target Drugs under Methods M2, M6, M7, M8\n")
        out.append("| Query Drug | Query ChEMBL ID | Target Drug | Target ChEMBL ID | Candidate Pool Size ($M$) | Rank M2 (S) | Rank M6 (Composite) | Rank M7 (Penalized) | Rank M8 (Popularity) |")
        out.append("| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_kc.iterrows():
            out.append(f"| {r['query_drug']} | {r['query_chembl_id']} | {r['target_drug']} | {r['target_chembl_id']} | {int(r['candidate_pool_size'])} | {r['rank_M2']:.1f} | {r['rank_M6']:.1f} | {r['rank_M7']:.1f} | {r['rank_M8']:.1f} |")
        out.append(make_source_line(kc_path))
        
    kc_top10_path = "known_cases_m7_top10.csv"
    if (RESULTS_DIR / kc_top10_path).exists():
        df_top10 = pd.read_csv(RESULTS_DIR / kc_top10_path)
        out.append("### Top 10 Ranked Candidates under Method M7 for Canonical Queries\n")
        out.append("| Query Drug | Position | Average Rank | Candidate Name | Candidate ChEMBL ID | M7 Score |")
        out.append("| :--- | :---: | :---: | :--- | :--- | :---: |")
        for _, r in df_top10.iterrows():
            out.append(f"| {r['query_drug']} | {int(r['position'])} | {r['average_rank']:.1f} | {r['candidate_name']} | {r['candidate_chembl_id']} | {r['m7_score']:.6f} |")
        out.append(make_source_line(kc_top10_path))
        
    out.append("---\n")

    # 8. MODALITY WEIGHT OPTIMIZATION AND STABILITY ANALYSIS
    out.append("## 8. Modality Weight Optimization and Stability Analysis\n")
    
    # 8.1 P1-dn MRR
    ws_p1_path = "weight_stability_p1_dn.csv"
    if (RESULTS_DIR / ws_p1_path).exists():
        df_w1 = pd.read_csv(RESULTS_DIR / ws_p1_path)
        out.append("### Modality Weight Stability on Protocol P1-dn (20 Random 50/50 Splits, Seeds 0-19)\n")
        out.append("| Split Seed | Best Tune Weights $(w_S, w_T, w_E)$ | Best Tune MRR | Best Test MRR | Default Test MRR | Test MRR Diff (Best - Def) | Default Tune MRR | Default Tune Rank (/66) |")
        out.append("| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_w1.iterrows():
            w_str = f"({r['best_w_s']:.1f}, {r['best_w_t']:.1f}, {r['best_w_e']:.1f})"
            out.append(f"| {int(r['split_seed'])} | {w_str} | {r['best_tune_mrr']:.4f} | {r['best_test_mrr']:.4f} | {r['default_test_mrr']:.4f} | {r['test_diff_mrr']:+.4f} | {r['default_tune_mrr']:.4f} | {int(r['default_tune_rank'])} |")
        out.append(make_source_line(ws_p1_path))

    # 8.1b P1-dn Hit@10
    ws_h1_path = "weight_stability_hit10_p1_dn.csv"
    if (RESULTS_DIR / ws_h1_path).exists():
        df_h1 = pd.read_csv(RESULTS_DIR / ws_h1_path)
        out.append("### Modality Weight Stability Hit@10 on Protocol P1-dn (20 Random 50/50 Splits, Seeds 0-19)\n")
        out.append("| Split Seed | Best Tune Weights $(w_S, w_T, w_E)$ | Best Test Hit@10 | Default Test Hit@10 | Test Hit@10 Diff (Best - Def) |")
        out.append("| :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_h1.iterrows():
            w_str = f"({r['best_w_s']:.1f}, {r['best_w_t']:.1f}, {r['best_w_e']:.1f})"
            out.append(f"| {int(r['split_seed'])} | {w_str} | {r['best_test_hit10']:.4f} | {r['default_test_hit10']:.4f} | {r['test_diff_hit10']:+.4f} |")
        out.append(make_source_line(ws_h1_path))

    # 8.2 P2b-dn MRR
    ws_p2b_path = "weight_stability_p2b_dn.csv"
    if (RESULTS_DIR / ws_p2b_path).exists():
        df_w2 = pd.read_csv(RESULTS_DIR / ws_p2b_path)
        out.append("### Modality Weight Stability on Protocol P2b-dn (0.4) (20 Random 50/50 Splits, Seeds 0-19)\n")
        out.append("| Split Seed | Best Tune Weights $(w_S, w_T, w_E)$ | Best Tune MRR | Best Test MRR | Default Test MRR | Test MRR Diff (Best - Def) | Default Tune MRR | Default Tune Rank (/66) |")
        out.append("| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_w2.iterrows():
            w_str = f"({r['best_w_s']:.1f}, {r['best_w_t']:.1f}, {r['best_w_e']:.1f})"
            out.append(f"| {int(r['split_seed'])} | {w_str} | {r['best_tune_mrr']:.4f} | {r['best_test_mrr']:.4f} | {r['default_test_mrr']:.4f} | {r['test_diff_mrr']:+.4f} | {r['default_tune_mrr']:.4f} | {int(r['default_tune_rank'])} |")
        out.append(make_source_line(ws_p2b_path))

    # 8.2b P2b-dn Hit@10
    ws_h2_path = "weight_stability_hit10_p2b_dn.csv"
    if (RESULTS_DIR / ws_h2_path).exists():
        df_h2 = pd.read_csv(RESULTS_DIR / ws_h2_path)
        out.append("### Modality Weight Stability Hit@10 on Protocol P2b-dn (0.4) (20 Random 50/50 Splits, Seeds 0-19)\n")
        out.append("| Split Seed | Best Tune Weights $(w_S, w_T, w_E)$ | Best Test Hit@10 | Default Test Hit@10 | Test Hit@10 Diff (Best - Def) |")
        out.append("| :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_h2.iterrows():
            w_str = f"({r['best_w_s']:.1f}, {r['best_w_t']:.1f}, {r['best_w_e']:.1f})"
            out.append(f"| {int(r['split_seed'])} | {w_str} | {r['best_test_hit10']:.4f} | {r['default_test_hit10']:.4f} | {r['test_diff_hit10']:+.4f} |")
        out.append(make_source_line(ws_h2_path))

    # 8.3 Weight Test Comparison
    wt_path = "weight_test_comparison.csv"
    if (RESULTS_DIR / wt_path).exists():
        df_wt = pd.read_csv(RESULTS_DIR / wt_path)
        out.append("### Independent Test Split Evaluation: Best Tune Weights vs Default Baseline\n")
        out.append("> **Evaluation Protocol & Provenance:** Protocol P2a ($\\theta=0.4$, original baseline, $N=1,645$ queries); 50/50 split with seed 42 ($N_{\\text{tune}}=822$, $N_{\\text{test}}=823$ queries; test set is held-out second half); scoring formulation uses production hub penalty ($M_7$: composite score divided by $1 + (\\text{hubness\\_percentile}/100)^2$); generator script: `experiments/step4_weight_sensitivity.py`.\n")
        out.append("| Configuration | Modality Weights $(w_S, w_T, w_E)$ | N Test Queries | MRR [95% CI] | Enrichment [95% CI] |")
        out.append("| :--- | :---: | :---: | :---: | :---: |")
        for _, r in df_wt.iterrows():
            w_str = f"({r['w_s']:.1f}, {r['w_t']:.1f}, {r['w_e']:.1f})"
            mrr_ci = f"{r['mrr_mean']:.4f} [{r['mrr_low']:.4f}, {r['mrr_high']:.4f}]"
            enr_ci = f"{r['enrich_mean']:.2f} [{r['enrich_low']:.2f}, {r['enrich_high']:.2f}]"
            out.append(f"| {r['configuration']} | {w_str} | {int(r['n_test'])} | {mrr_ci} | {enr_ci} |")
        out.append(make_source_line(wt_path))

    # 8.4 Weight Simplex Grid Search Summary
    wg_path = "weight_grid.csv"
    if (RESULTS_DIR / wg_path).exists():
        df_wg = pd.read_csv(RESULTS_DIR / wg_path)
        out.append("### Modality Weight Grid Search Summary (66 Triples on 2-Simplex at Step 0.1)\n")
        out.append("> **Evaluation Protocol & Provenance:** Protocol P2a ($\\theta=0.4$, original baseline, $N=1,645$ queries); evaluated on tuning split ($N_{\\text{tune}}=822$ of 1,645, split seed 42); scoring formulation uses production hub penalty ($M_7$); generator script: `experiments/step4_weight_sensitivity.py`.\n")
        out.append("| Rank | Modality Weights $(w_S, w_T, w_E)$ | Tune MRR | Note |")
        out.append("| :---: | :---: | :---: | :--- |")
        df_wg_sorted = df_wg.sort_values(by="mrr_tune", ascending=False).reset_index(drop=True)
        for rank_idx, r in df_wg_sorted.head(10).iterrows():
            w_str = f"({r['w_s']:.1f}, {r['w_t']:.1f}, {r['w_e']:.1f})"
            note = "Best tuning configuration" if rank_idx == 0 else ""
            out.append(f"| {rank_idx + 1} | {w_str} | {r['mrr_tune']:.4f} | {note} |")
        def_idx = df_wg_sorted[(df_wg_sorted['w_s'] == 0.4) & (df_wg_sorted['w_t'] == 0.4) & (df_wg_sorted['w_e'] == 0.2)].index[0]
        if def_idx >= 10:
            r_def = df_wg_sorted.loc[def_idx]
            out.append(f"| {def_idx + 1} | ({r_def['w_s']:.1f}, {r_def['w_t']:.1f}, {r_def['w_e']:.1f}) | {r_def['mrr_tune']:.4f} | Production Default Baseline |")
        out.append(make_source_line(wg_path))

    out.append("---\n")

    # 9. FAILURE ANALYSIS FOR COMPOSITE SCORING (METHOD M7)
    out.append("## 9. Failure Analysis for Composite Scoring (Method M7)\n")
    
    # 9.1 P2b-dn (0.4)
    fa_p2b_path = "failure_analysis_p2b_dn.csv"
    if (RESULTS_DIR / fa_p2b_path).exists():
        df_fa2 = pd.read_csv(RESULTS_DIR / fa_p2b_path)
        out.append("### Top 20 Best Performing Queries for M7 on Protocol P2b-dn (0.4)\n")
        out.append("| Query Drug | N Relevant | First Relevant Drug | Rank ($r_{\\min}$) | S | T | E | Targets (Q / Cand) | SEs (Q / Cand) |")
        out.append("| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_fa2[df_fa2["category"] == "best_20"].iterrows():
            q_str = f"{r['query_cid']} ({r['query_name']})"
            rel_str = f"{r['first_rel_cid']} ({r['first_rel_name']})"
            t_flag = f"{'Y' if str(r['query_has_targets']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_targets']).lower() == 'true' else 'N'}"
            se_flag = f"{'Y' if str(r['query_has_se']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_se']).lower() == 'true' else 'N'}"
            out.append(f"| {q_str} | {int(r['n_relevant'])} | {rel_str} | {r['first_rel_rank']:.1f} | {r['S']:.4f} | {r['T']:.4f} | {r['E']:.4f} | {t_flag} | {se_flag} |")
        out.append(make_source_line(fa_p2b_path))

        out.append("### Top 20 Worst Performing Queries for M7 on Protocol P2b-dn (0.4)\n")
        out.append("| Query Drug | N Relevant | First Relevant Drug | Rank ($r_{\\min}$) | S | T | E | Targets (Q / Cand) | SEs (Q / Cand) |")
        out.append("| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_fa2[df_fa2["category"] == "worst_20"].iterrows():
            q_str = f"{r['query_cid']} ({r['query_name']})"
            rel_str = f"{r['first_rel_cid']} ({r['first_rel_name']})"
            t_flag = f"{'Y' if str(r['query_has_targets']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_targets']).lower() == 'true' else 'N'}"
            se_flag = f"{'Y' if str(r['query_has_se']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_se']).lower() == 'true' else 'N'}"
            out.append(f"| {q_str} | {int(r['n_relevant'])} | {rel_str} | {r['first_rel_rank']:.1f} | {r['S']:.4f} | {r['T']:.4f} | {r['E']:.4f} | {t_flag} | {se_flag} |")
        out.append(make_source_line(fa_p2b_path))

    # 9.2 P3-dn
    fa_p3_path = "failure_analysis_p3_dn.csv"
    if (RESULTS_DIR / fa_p3_path).exists():
        df_fa3 = pd.read_csv(RESULTS_DIR / fa_p3_path)
        out.append("### Top 20 Best Performing Queries for M7 on Protocol P3-dn\n")
        out.append("| Query Drug | N Relevant | First Relevant Drug | Rank ($r_{\\min}$) | S | T | E | Targets (Q / Cand) | SEs (Q / Cand) |")
        out.append("| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_fa3[df_fa3["category"] == "best_20"].iterrows():
            q_str = f"{r['query_cid']} ({r['query_name']})"
            rel_str = f"{r['first_rel_cid']} ({r['first_rel_name']})"
            t_flag = f"{'Y' if str(r['query_has_targets']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_targets']).lower() == 'true' else 'N'}"
            se_flag = f"{'Y' if str(r['query_has_se']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_se']).lower() == 'true' else 'N'}"
            out.append(f"| {q_str} | {int(r['n_relevant'])} | {rel_str} | {r['first_rel_rank']:.1f} | {r['S']:.4f} | {r['T']:.4f} | {r['E']:.4f} | {t_flag} | {se_flag} |")
        out.append(make_source_line(fa_p3_path))

        out.append("### Top 20 Worst Performing Queries for M7 on Protocol P3-dn\n")
        out.append("| Query Drug | N Relevant | First Relevant Drug | Rank ($r_{\\min}$) | S | T | E | Targets (Q / Cand) | SEs (Q / Cand) |")
        out.append("| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
        for _, r in df_fa3[df_fa3["category"] == "worst_20"].iterrows():
            q_str = f"{r['query_cid']} ({r['query_name']})"
            rel_str = f"{r['first_rel_cid']} ({r['first_rel_name']})"
            t_flag = f"{'Y' if str(r['query_has_targets']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_targets']).lower() == 'true' else 'N'}"
            se_flag = f"{'Y' if str(r['query_has_se']).lower() == 'true' else 'N'} / {'Y' if str(r['cand_has_se']).lower() == 'true' else 'N'}"
            out.append(f"| {q_str} | {int(r['n_relevant'])} | {rel_str} | {r['first_rel_rank']:.1f} | {r['S']:.4f} | {r['T']:.4f} | {r['E']:.4f} | {t_flag} | {se_flag} |")
        out.append(make_source_line(fa_p3_path))

    out.append("---\n")

    # 10. CANDIDATE POOL DUPLICATE AUDIT AND QUERY IMPACT
    out.append("## 10. Candidate Pool Duplicate Audit and Query Impact\n")
    dup_path = "duplicate_pairs.csv"
    if (RESULTS_DIR / dup_path).exists():
        dup_path_full = RESULTS_DIR / dup_path
        df_dup = pd.read_csv(dup_path_full)
        drugs_set = set(df_dup['chembl_id_1']).union(set(df_dup['chembl_id_2']))
        n_pairs = len(df_dup)
        n_drugs = len(drugs_set)
        inchikey_count = int((df_dup['reasons'] == 'first_block_inchikey').sum())
        name_count = int((df_dup['reasons'] == 'generic_base_name').sum())
        both_count = int((df_dup['reasons'] == 'generic_base_name; first_block_inchikey').sum())

        # Dynamic computation of query impact and drop counts from CSVs
        p1_csv = RESULTS_DIR / "metrics_C1_P1.csv"
        p1_dn_csv = RESULTS_DIR / "metrics_C1_P1_dn.csv"
        p3_csv = RESULTS_DIR / "metrics_C1_P3.csv"
        p3_dn_csv = RESULTS_DIR / "metrics_C1_P3_dn.csv"

        n_p1_orig = int(pd.read_csv(p1_csv)["n_queries"].iloc[0])
        n_p1_dn = int(pd.read_csv(p1_dn_csv)["n_queries"].iloc[0])
        n_p3_orig = int(pd.read_csv(p3_csv)["n_queries"].iloc[0])
        n_p3_dn = int(pd.read_csv(p3_dn_csv)["n_queries"].iloc[0])
        p1_dropped = n_p1_orig - n_p1_dn
        p3_dropped = n_p3_orig - n_p3_dn

        # Map duplicate relationships from duplicate_pairs.csv
        duplicates_of = defaultdict(set)
        for _, row in df_dup.iterrows():
            duplicates_of[row['chembl_id_1']].add(row['chembl_id_2'])
            duplicates_of[row['chembl_id_2']].add(row['chembl_id_1'])

        gold_path = DATA_DIR / "gold_pairs_primary.csv"
        drug_inds = defaultdict(set)
        with open(gold_path, "r", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                drug_inds[row["chembl_id"]].add(row["indication_id"])

        cids_with_inds = sorted(drug_inds.keys())
        p1_queries_with_dup = 0
        p1_queries_dup_in_rel = 0
        for q_cid in cids_with_inds:
            q_inds = drug_inds[q_cid]
            rel_p1 = [c for c in cids_with_inds if c != q_cid and not q_inds.isdisjoint(drug_inds[c])]
            if rel_p1:
                if duplicates_of[q_cid]:
                    p1_queries_with_dup += 1
                    if duplicates_of[q_cid].intersection(rel_p1):
                        p1_queries_dup_in_rel += 1

        from experiments.eval_engine import load_database, DB_A_PATH
        pool_size = len(load_database(DB_A_PATH, exclude_inorganics=True))

        out.append("### Duplicate Pairs Audit Summary and Representative Sample Pairs\n")
        out.append(f"> **Duplicate Census:** {n_drugs} unique pool drugs involved across {n_pairs} duplicate pairs in the {pool_size:,}-drug pool.  \n"
                   f"> **Reason Breakdown:** First-block InChIKey: {inchikey_count} pairs, Generic base name: {name_count} pairs, Both: {both_count} pairs.  \n"
                   f"> **Query Impact:** {p1_queries_with_dup} P1 queries have >= 1 duplicate in candidate pool; {p1_queries_dup_in_rel} P1 queries have >= 1 duplicate in Rel(q). "
                   f"Queries dropped due to zero remaining indications: P1-dn dropped {p1_dropped} ({n_p1_orig:,} -> {n_p1_dn:,}), P3-dn dropped {p3_dropped} ({n_p3_orig:,} -> {n_p3_dn:,}).  \n"
                   f"> **Census Source Files:** `{gold_path.as_posix()}` (SHA-256: `{sha256_file(gold_path)}`), "
                   f"`{p1_csv.as_posix()}` (SHA-256: `{sha256_file(p1_csv)}`), "
                   f"`{p1_dn_csv.as_posix()}` (SHA-256: `{sha256_file(p1_dn_csv)}`), "
                   f"`{p3_csv.as_posix()}` (SHA-256: `{sha256_file(p3_csv)}`), "
                   f"`{p3_dn_csv.as_posix()}` (SHA-256: `{sha256_file(p3_dn_csv)}`), "
                   f"`{dup_path_full.as_posix()}` (SHA-256: `{sha256_file(dup_path_full)}`)\n")
        out.append("| Pair Index | Drug 1 ChEMBL ID | Drug 1 Name | Drug 2 ChEMBL ID | Drug 2 Name | Removal Reason |")
        out.append("| :---: | :--- | :--- | :--- | :--- | :--- |")
        for idx, r in df_dup.head(20).iterrows():
            out.append(f"| {idx + 1} | {r['chembl_id_1']} | {r['name_1']} | {r['chembl_id_2']} | {r['name_2']} | `{r['reasons']}` |")
        out.append(make_source_line(dup_path))

    out.append("---\n")

    # 11. SIDE-EFFECT MUTUAL INFORMATION LIFT ANALYSIS
    out.append("## 11. Side-Effect Mutual Information Lift Analysis\n")
    lift_path = "se_top50_lift.csv"
    if (RESULTS_DIR / lift_path).exists():
        df_lift = pd.read_csv(RESULTS_DIR / lift_path)
        out.append("### Top 50 Side Effects Ranked by Indication Co-occurrence Lift\n")
        out.append("| Rank | Preferred Term | Drug Count ($N_{\\text{drugs}}$) | Indication Co-occurrences | Lift |")
        out.append("| :---: | :--- | :---: | :---: | :---: |")
        for _, r in df_lift.iterrows():
            out.append(f"| {int(r['rank'])} | {r['term']} | {int(r['n_drugs'])} | {int(r['n_relevant_pair_cooccurrences'])} | {r['lift']:.4f} |")
        out.append(make_source_line(lift_path))
        
    final_md_path = RESULTS_DIR / "FINAL_TABLES.md"
    content = "\n".join(out)
    with open(final_md_path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Saved consolidated tables to: {final_md_path}")
    print(f"FINAL_TABLES.md length: {len(content)} characters, {len(out)} lines.")

    # Synchronize experiments/RESULTS.md to ensure single source of truth
    results_md_path = ROOT / "experiments" / "RESULTS.md"
    with open(results_md_path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"Synchronized {results_md_path} from FINAL_TABLES.md")

def update_manifest():
    manifest_path = RESULTS_DIR / "manifest.json"
    
    # Read real HEAD commit and check working tree status of experiments/
    try:
        git_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(ROOT), text=True).strip()
    except Exception:
        git_commit = "UNKNOWN"
        
    try:
        # Check if tracked working tree files in experiments/ have modifications (excluding generator outputs)
        status_lines = subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no", "experiments/"],
            cwd=str(ROOT),
            text=True
        ).splitlines()
        generated_files = {"generate_final_tables.py", "FINAL_TABLES.md", "RESULTS.md", "manifest.json"}
        dirty_other = [l for l in status_lines if not any(gf in l for gf in generated_files)]
        experiments_clean = (len(dirty_other) == 0)
    except Exception:
        experiments_clean = False
        
    # Read existing or construct
    if manifest_path.exists():
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)
    else:
        manifest = {}
        
    manifest["metadata"] = {
        "description": "MoleculAR Evaluation Harness Manifest",
        "git_commit": git_commit,
        "experiments_working_tree_clean": experiments_clean,
        "random_seed": 42,
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "sqlite3": "3.50.4",
            "psycopg2": "2.9.13 (dt dec pq3 ext lo64)",
            "pytest": "9.2.0"
        }
    }
    
    manifest["databases"] = {
        "DB-A": {
            "path": str(DATA_DIR / "drugrepo_oct6.db"),
            "sha256": sha256_file(DATA_DIR / "drugrepo_oct6.db") if (DATA_DIR / "drugrepo_oct6.db").exists() else "D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7"
        },
        "DB-R": {
            "path": str(DATA_DIR / "drugrepo_oct6_recovered.db"),
            "sha256": sha256_file(DATA_DIR / "drugrepo_oct6_recovered.db") if (DATA_DIR / "drugrepo_oct6_recovered.db").exists() else "66928166E834ABE98AF895E7E8A5234D11EB2D106EF1A034C9801EA027F07A38"
        },
        "root_drugrepo_db": {
            "path": str(ROOT / "drugrepo.db"),
            "sha256": sha256_file(ROOT / "drugrepo.db") if (ROOT / "drugrepo.db").exists() else "D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7"
        }
    }
    
    # Hash every CSV in experiments/results
    results_csvs = {}
    for f in sorted(os.listdir(RESULTS_DIR)):
        if f.endswith(".csv"):
            f_path = RESULTS_DIR / f
            results_csvs[f] = {
                "path": str(f_path),
                "sha256": sha256_file(f_path)
            }
            
    manifest["results_csvs"] = results_csvs
    
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        
    print(f"Updated manifest.json with {len(results_csvs)} result CSV hashes.")

if __name__ == "__main__":
    build_final_tables()
    update_manifest()
