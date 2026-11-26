"""
experiments/generate_manifest.py

Generates experiments/results/manifest.json containing:
- SHA-256 of DB-A, DB-R, root DB, all gold files, rejected files, and every results CSV.
- Git commit hash (read-only git rev-parse HEAD).
- Random seed (42).
- Library and Python versions.
- Exact CLI commands to reproduce every step.
"""

import hashlib
import json
import subprocess
import sys
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "experiments" / "data"
RESULTS_DIR = ROOT / "experiments" / "results"
MANIFEST_PATH = RESULTS_DIR / "manifest.json"

def get_sha256(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest().upper()

def get_git_commit() -> str:
    res = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=ROOT)
    return res.stdout.strip()

def main():
    import numpy as np
    import rdkit
    import psycopg2
    import matplotlib
    import pytest

    # Databases
    databases = {
        "DB-A": {
            "path": str(DATA_DIR / "drugrepo_oct6.db"),
            "sha256": get_sha256(DATA_DIR / "drugrepo_oct6.db"),
        },
        "DB-R": {
            "path": str(DATA_DIR / "drugrepo_oct6_recovered.db"),
            "sha256": get_sha256(DATA_DIR / "drugrepo_oct6_recovered.db"),
        },
        "root_drugrepo_db": {
            "path": str(ROOT / "drugrepo.db"),
            "sha256": get_sha256(ROOT / "drugrepo.db"),
        },
    }

    # Gold standards & input data
    gold_files = {}
    for p in [
        DATA_DIR / "gold_pairs_primary.csv",
        DATA_DIR / "gold_pairs_rule_a_only.csv",
        DATA_DIR / "gold_pairs_with_inorganics.csv",
        DATA_DIR / "excluded_indications.csv",
        DATA_DIR / "rejected_rule_b.csv",
        DATA_DIR / "drugcentral_raw_atc.csv",
        DATA_DIR / "drugcentral_raw_approved_indications.csv",
    ]:
        if p.exists():
            gold_files[p.name] = {
                "path": str(p),
                "sha256": get_sha256(p),
            }

    # Result CSVs
    results_csvs = {}
    for p in sorted(RESULTS_DIR.glob("*.csv")):
        results_csvs[p.name] = {
            "path": str(p),
            "sha256": get_sha256(p),
        }

    manifest = {
        "metadata": {
            "description": "MoleculAR Evaluation Harness Manifest",
            "git_commit": get_git_commit(),
            "code_commit": get_git_commit(),
            "random_seed": 42,
            "environment": {
                "python": sys.version,
                "numpy": np.__version__,
                "rdkit": rdkit.__version__,
                "sqlite3": sqlite3.sqlite_version,
                "psycopg2": psycopg2.__version__,
                "matplotlib": matplotlib.__version__,
                "pytest": pytest.__version__,
            },
        },
        "databases": databases,
        "gold_and_reference_files": gold_files,
        "results_csvs": results_csvs,
        "reproduction_commands": {
            "step1_rule_b_and_gold": "python experiments/step1_update_rule_b.py",
            "step2_reproduce_production": "python experiments/step2_reproduce_production.py",
            "step3_run_all_methods": "python experiments/step3_run_all_methods.py",
            "step4_weight_sensitivity": "python experiments/step4_weight_sensitivity.py",
            "step5_pytest": "pytest experiments/tests/test_harness.py",
            "step5_verify_experiments": "python experiments/verify_experiments.py",
            "step6_visualizations": "python experiments/step6_visualizations.py",
            "step6_generate_results_md": "python experiments/generate_results_md.py",
            "step6_generate_manifest": "python experiments/generate_manifest.py",
            "run_all": "python experiments/run_all.py",
        },
    }

    with open(MANIFEST_PATH, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"Generated manifest.json -> {MANIFEST_PATH}")

if __name__ == "__main__":
    main()

