# REPRODUCIBILITY GUIDE: MoleculAR Evaluation Harness

This document provides exact, sequential instructions to regenerate every experimental result, benchmark metric, null model, ablation, and finding from the raw primary database **DB-A** (`drugrepo_oct6.db`).

---

## 1. System Requirements & Environment

- **Operating System:** Windows, Linux, or macOS
- **Python Version:** Python 3.10+ (Tested on Python 3.14.5)
- **Key Dependencies:**
  ```bash
  pip install numpy scipy pandas rdkit psycopg2-binary pytest
  ```
- **Random Seed:** `42` (Fixed across all bootstrap resamplings, permutations, and sampling steps)
- **Floating-Point Precision:** IEEE-754 64-bit float (`np.float64`)
- **Tie-Breaking Governance:** Exact mid-rank average ranking:
  $$r(c) = 1 + N_{>} + 0.5 \times N_{=}$$

---

## 2. Primary Data Integrity Check

Before running any experiment, verify that the primary database is present and uncorrupted:

```bash
# Verify DB-A file integrity and page-level invariants
python -c "
import hashlib
p = 'experiments/data/drugrepo_oct6.db'
sha = hashlib.sha256(open(p, 'rb').read()).hexdigest().upper()
print('DB-A SHA-256:', sha)
assert sha == 'D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7', 'Hash mismatch!'
print('Integrity check PASSED.')
"
```

---

## 3. Step-by-Step Regeneration Pipeline

Execute the following commands in exact numerical sequence from the repository root:

### Step 1: Gold Standard Indication Pairs
Extract therapeutic indications from DB-A, map to ChEMBL IDs, apply stereo-prefix exclusions to Rule (b), and output gold pair files:
```bash
python experiments/step1_gold_standard.py
```
- **Outputs:**
  - `experiments/data/gold_pairs_primary.csv`
  - `experiments/data/gold_pairs_rule_a_only.csv`
  - `experiments/data/rejected_rule_b.csv`

---

### Step 2: Primary Evaluation Engine (Protocols P1, P2a, P2b, P3)
Precompute pairwise similarity matrices ($S, T, E$), hubness penalties, and evaluate Methods M1 through M9 across Candidate Sets C1–C4:
```bash
python experiments/step3_run_all_methods.py
```
- **Outputs:**
  - `experiments/results/metrics_C{1,2,3,4}_P{1,2a_0.3,2a_0.4,2a_0.5,2b_0.4,3,3_all}.csv`

---

### Step 3: Duplicate Audit and Duplicate-Neutral Protocols (-dn)
Audit same-drug duplicates (generic base name, first-block InChIKey, or desalted canonical SMILES), run duplicate-neutral protocols, and paired comparisons:
```bash
python experiments/step1_duplicate_neutral_run.py
python experiments/step2_paired_comparisons.py
python experiments/step5_failure_analysis.py
```
- **Outputs:**
  - `experiments/results/duplicate_pairs.csv`
  - `experiments/results/metrics_C1_P1_dn.csv`
  - `experiments/results/metrics_C1_P2a_dn_0.4.csv`
  - `experiments/results/metrics_C1_P2b_dn_0.4.csv`
  - `experiments/results/metrics_C1_P3_dn.csv`
  - `experiments/results/paired_differences.csv`
  - `experiments/results/failure_analysis_p2b_dn.csv`
  - `experiments/results/failure_analysis_p3_dn.csv`

---

### Step 4: Coverage-Stratified Evaluation (Protocols P1, P2a, P2b, P3-dn)
Stratify queries by data availability (Targets + SE, Targets only, SE only, Neither) and evaluate Methods M2–M7:
```bash
python experiments/step7_coverage_stratification.py
```
- **Outputs:**
  - `experiments/results/coverage_stratified_all.csv`
  - `experiments/results/coverage_stratified_P1.csv`
  - `experiments/results/coverage_stratified_P2a_0.4.csv`
  - `experiments/results/coverage_stratified_P2b_0.4.csv`
  - `experiments/results/coverage_stratified_P3_dn.csv`

---

### Step 5: Stricter Cross-Class Protocols P4 and P4b
Extract WHO ATC Level-1 anatomical main group letters, audit cross-organ indications, and run duplicate-neutral Protocols P4 and P4b:
```bash
python experiments/step11_p4_protocol.py
python experiments/step14_run_p4_p4b.py
```
- **Outputs:**
  - `experiments/results/p4_query_census.csv`
  - `experiments/results/p4_p4b_query_census.csv`
  - `experiments/results/metrics_C1_P4.csv`
  - `experiments/results/metrics_C1_P4b.csv`
  - `experiments/results/paired_differences_p4_p4b.csv`
  - `experiments/results/coverage_stratified_p4_p4b.csv`

---

### Step 6: Null Models and Feature Shuffle Controls
Run 20 feature-shuffled controls (seed 42) to establish empirical null distributions:
```bash
# Unstratified shuffle controls
python experiments/step10_shuffle_control.py

# Coverage-fair within-strata shuffle controls
python experiments/step15_coverage_fair_nulls.py
```
- **Outputs:**
  - `experiments/results/shuffle_control_results.csv`
  - `experiments/results/coverage_fair_nulls.csv`

---

### Step 7: Ablation Studies & Sensitivity Analyses
```bash
# Hubness penalty sensitivity sweep (k in [0.0, 3.0])
python experiments/step4_hub_ablation.py

# Weight grid search (66 triples) and test comparison:
python experiments/step4_weight_sensitivity.py

# Weight stability analysis (20 random splits)
python experiments/step3_weight_stability.py

# Baseline side effect term co-occurrence and MI leakage analysis:
python experiments/step8_se_leakage.py

# Side effect top-50 lift audit
python experiments/step12_se_lift_terms.py

# Approved side-effect leakage ablations (Conditions A, B, C)
python experiments/step16_se_leakage_ablations.py

# Tie sensitivity analysis (Average vs Best vs Worst rank)
python experiments/step18_tie_sensitivity.py
python experiments/step20_m8_tie_sensitivity.py

# Canonical 8 known-case checks on DB-A
python experiments/step19_known_cases_db_a.py
```
- **Outputs:**
  - `experiments/results/hub_penalty_ablation.csv`
  - `experiments/results/weight_grid.csv`
  - `experiments/results/weight_test_comparison.csv`
  - `experiments/results/weight_stability_p1_dn.csv`
  - `experiments/results/weight_stability_p2b_dn.csv`
  - `experiments/results/se_leakage_terms.csv` (Note: `se_leakage_terms.csv` may differ in the order of two exactly tied terms, ranks 20/21)
  - `experiments/results/se_leakage_ablation.csv`
  - `experiments/results/se_top50_lift.csv`
  - `experiments/results/se_leakage_approved_ablations.csv`
  - `experiments/results/tie_sensitivity.csv`
  - `experiments/results/tie_sensitivity_m8.csv`
  - `experiments/results/known_cases_db_a.csv`
  - `experiments/results/known_cases_m7_top10.csv`

---

## 4. Verification and Test Suite

Execute the independent brute-force verifiers and pytest suites:

```bash
# 1. Independent brute-force rank verification for M2-M9 (Step 5.3)
python experiments/verify_experiments.py

# 2. Independent brute-force rank verification for Protocols P4 and P4b
python experiments/verify_p4_p4b.py

# 3. Main harness unit test suite
pytest experiments/tests/test_harness.py -v

# 4. P4/P4b brute-force and P1 exact reproduction pytest suite
pytest experiments/verify_p4_p4b.py -v
```

---

## 5. Report & Table Consolidation

Finally, compile all CSV results into unified markdown tables and update the manifest:

```bash
# Generate second-round findings tables
python experiments/generate_findings_tables.py

# Generate final consolidated tables and update manifest.json
python experiments/generate_final_tables.py
```
- **Outputs:**
  - `experiments/results/FINDINGS_TABLES.md`
  - `experiments/results/FINAL_TABLES.md`
  - `experiments/results/manifest.json`
