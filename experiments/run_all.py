"""
experiments/run_all.py

Master Orchestrator for MoleculAR Evaluation Harness:
Regenerates the complete evaluation suite from scratch:
1. Step 1: Rule (b) stereo-prefix rejection, gold standard regeneration.
2. Step 2: Production reproduction verification against restored pipeline.py and root drugrepo.db.
3. Build C3 and C4 auxiliary gold standards.
4. Step 3: Run all methods (M1-M9) across all configurations (C1-C4) and protocols, plus C1 size strata.
5. Step 4: Weight sensitivity grid search (66 triples) and heatmap.
6. Step 5: Unit test suite (pytest) and brute-force ranking verification (verify_experiments.py).
7. Step 6: Visualizations (MRR and Enrichment bar charts).
8. Step 6: Script-generated RESULTS.md.
9. Step 6: manifest.json with all SHAs, commit hash, seed, and library versions.
"""

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

def run_step(cmd_list, description):
    print(f"\n{'='*80}")
    print(f"RUNNING: {description}")
    print(f"COMMAND: {' '.join(cmd_list)}")
    print(f"{'='*80}")
    t0 = time.time()
    res = subprocess.run(cmd_list, cwd=ROOT)
    if res.returncode != 0:
        print(f"ERROR: Step failed with return code {res.returncode}: {description}")
        sys.exit(res.returncode)
    print(f"COMPLETED in {time.time() - t0:.2f}s: {description}\n")

def main():
    start_time = time.time()
    py = sys.executable

    print("\n" + "#"*80)
    print("STARTING MOLECULAR EVALUATION HARNESS: END-TO-END REGENERATION")
    print("#"*80 + "\n")

    # Step 1: Rule (b) update & Gold regeneration
    run_step([py, "experiments/step1_update_rule_b.py"], "Step 1: Rule (b) Rejections & Gold Standard")

    # Step 2: Production reproduction check
    run_step([py, "experiments/step2_reproduce_production.py"], "Step 2: Production Reproduction Check")

    # C3 & C4 Gold files
    run_step([py, "experiments/build_c3_c4_gold.py"], "Build C3 & C4 Gold Standards")

    # Step 3: Run all methods across all configs and protocols
    run_step([py, "experiments/step3_run_all_methods.py"], "Step 3: Run All Methods (C1-C4, P1-P3, Size Strata)")

    # Step 4: Weight sensitivity grid search
    run_step([py, "experiments/step4_weight_sensitivity.py"], "Step 4: Weight Sensitivity Grid Search & Heatmap")

    # Step 5: Pytest unit tests
    run_step([py, "-m", "pytest", "experiments/tests/test_harness.py", "-v"], "Step 5.1 & 5.2: Pytest Suite & Toy Database")

    # Step 5: Brute-force verification
    run_step([py, "experiments/verify_experiments.py"], "Step 5.3: Independent Brute-Force Ranking Verification")

    # Step 6: Visualizations
    run_step([py, "experiments/step6_visualizations.py"], "Step 6: Generate MRR & Enrichment Visualizations")

    # Step 6: RESULTS.md
    run_step([py, "experiments/generate_results_md.py"], "Step 6: Generate RESULTS.md from CSVs")

    # Step 6: manifest.json
    run_step([py, "experiments/generate_manifest.py"], "Step 6: Generate manifest.json")

    print("\n" + "#"*80)
    print(f"ALL STEPS COMPLETED SUCCESSFULLY IN {time.time() - start_time:.2f}s!")
    print("#"*80 + "\n")

if __name__ == "__main__":
    main()

