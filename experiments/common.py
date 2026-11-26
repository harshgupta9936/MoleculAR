"""
experiments/common.py — Shared data loader for the MoleculAR evaluation harness.

Loads the snapshot database and builds per-chembl_id feature dictionaries,
matching the production code's table reads (drugs, targets, side_effects).
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path

from rdkit.Chem import DataStructs


def load_database(db_path: str) -> dict:
    """
    Load the drugrepo snapshot database and return a dict keyed by chembl_id.

    Each value is a dict with keys:
        chembl_id, generic_base_name, drug_name, smiles, indication,
        fingerprint (rdkit ExplicitBitVect or None),
        targets (set of uniprot_id strings),
        side_effects (set of side_effect_name strings),
        hubness_percentile (float or None),
        is_trivial (bool),
        trivial_category (str or None)
    """
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    drugs = {}

    # ---- Load drugs table (same columns production code reads) ----
    cur.execute(
        "SELECT chembl_id, generic_base_name, drug_name, smiles, indication, "
        "fingerprint_blob, hubness_percentile, trivial_category FROM drugs"
    )
    for row in cur.fetchall():
        cid, gbn, dname, smiles, indication, fp_blob, hub_pct, trivial = row
        fp = DataStructs.CreateFromBitString(fp_blob) if fp_blob else None
        drugs[cid] = {
            "chembl_id": cid,
            "generic_base_name": gbn,
            "drug_name": dname,
            "smiles": smiles,
            "indication": indication,
            "fingerprint": fp,
            "targets": set(),
            "side_effects": set(),
            "hubness_percentile": hub_pct,
            "is_trivial": trivial is not None,
            "trivial_category": trivial,
        }

    # ---- Load targets (same join as production) ----
    cur.execute("SELECT chembl_id, uniprot_id FROM targets")
    for cid, uid in cur.fetchall():
        if cid in drugs:
            drugs[cid]["targets"].add(uid)

    # ---- Load side effects (same join as production) ----
    cur.execute("SELECT chembl_id, side_effect_name FROM side_effects")
    for cid, se in cur.fetchall():
        if cid in drugs:
            drugs[cid]["side_effects"].add(se)

    conn.close()
    return drugs


def scored_drugs(drugs: dict) -> dict:
    """Return only non-trivial drugs (the scored set)."""
    return {cid: d for cid, d in drugs.items() if not d["is_trivial"]}


def print_census(drugs: dict):
    """Print and return census counts."""
    all_ids = set(drugs.keys())
    nontrivial = {cid for cid, d in drugs.items() if not d["is_trivial"]}
    with_targets = {cid for cid in nontrivial if len(drugs[cid]["targets"]) >= 1}
    with_se = {cid for cid in nontrivial if len(drugs[cid]["side_effects"]) >= 1}
    with_both = with_targets & with_se
    with_neither = nontrivial - with_targets - with_se

    counts = {
        "total_drugs": len(all_ids),
        "non_trivial": len(nontrivial),
        "with_targets": len(with_targets),
        "with_side_effects": len(with_se),
        "with_both": len(with_both),
        "with_neither": len(with_neither),
    }

    for k, v in counts.items():
        print(f"  {k}: {v}")

    return counts


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="Path to drugrepo_snapshot.db")
    return p.parse_args()


if __name__ == "__main__":
    args = get_args()
    print(f"Loading database from {args.db} ...")
    drugs = load_database(args.db)
    print("Census:")
    counts = print_census(drugs)

    # Save census to logs
    log_dir = Path(args.db).parent.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    census_path = log_dir / "census.json"
    with open(census_path, "w") as f:
        json.dump(counts, f, indent=2)
    print(f"Saved census to {census_path}")

