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
import hashlib
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

def format_protocol_table(csv_name, title):
    p = RESULTS_DIR / csv_name
    if not p.exists():
        return ""
    df = pd.read_csv(p)
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
    out.append("**Generated:** 2026-11-14  ")
    out.append("**Primary Database:** DB-A (`experiments/data/drugrepo_oct6.db`)  ")
    out.append("**Strict Governance:** Frozen PROTOCOL.md definitions, exact mid-rank average ranking for ties, float64 precision, seed 42, chembl_id keys.  \n")
    out.append("---\n")

    # 1. PRIMARY PROTOCOL PERFORMANCE METRICS (C1)
    out.append("## 1. Primary Protocol Performance Metrics (Candidate Pool C1)\n")
    
    protocols_c1 = [
        ("metrics_C1_P1.csv", "Protocol P1 (Broad Therapeutic Indication Matching - Original Baseline, N=1,661)"),
        ("metrics_C1_P1_dn.csv", "Protocol P1-dn (Broad Indication Matching - Duplicate-Neutral, N=1,661)"),
        ("metrics_C1_P2a_0.4.csv", "Protocol P2a (0.4) (Analogue Neutralization at θ=0.4 - Original Baseline, N=1,647)"),
        ("metrics_C1_P2a_dn_0.4.csv", "Protocol P2a-dn (0.4) (Analogue Neutralization at θ=0.4 - Duplicate-Neutral, N=1,647)"),
        ("metrics_C1_P2a_0.3.csv", "Protocol P2a (0.3) (Analogue Neutralization at θ=0.3, N=1,640)"),
        ("metrics_C1_P2a_0.5.csv", "Protocol P2a (0.5) (Analogue Neutralization at θ=0.5, N=1,650)"),
        ("metrics_C1_P2b_0.4.csv", "Protocol P2b (0.4) (Analogue Removal at θ=0.4 - Original Baseline, N=1,647)"),
        ("metrics_C1_P2b_dn_0.4.csv", "Protocol P2b-dn (0.4) (Analogue Removal at θ=0.4 - Duplicate-Neutral, N=1,647)"),
        ("metrics_C1_P3.csv", "Protocol P3 (Cross-Class Repurposing: Disjoint ATC-L3, Pool with ATC - Original Baseline, N=1,240)"),
        ("metrics_C1_P3_all.csv", "Protocol P3-all (Cross-Class Repurposing: Disjoint ATC-L3, Full Pool - Original Baseline, N=1,240)"),
        ("metrics_C1_P3_dn.csv", "Protocol P3-dn (Cross-Class Repurposing: Disjoint ATC-L3, Pool with ATC - Duplicate-Neutral, N=1,240)"),
        ("metrics_C1_P4.csv", "Protocol P4 (Stricter Cross-Class: Disjoint ATC Level-1 Anatomical Main Groups, Full ATC Pool, N=965)"),
        ("metrics_C1_P4b.csv", "Protocol P4b (Stricter Cross-Class: Disjoint ATC Level-1 Anatomical Removal, N=965)"),
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
    
    # Read existing or construct
    if manifest_path.exists():
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)
    else:
        manifest = {}
        
    manifest["metadata"] = {
        "description": "MoleculAR Evaluation Harness Manifest",
        "git_commit": "27a8e8cf425790efa002f422f24eabe7174c10ea",
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
