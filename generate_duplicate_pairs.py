"""
generate_duplicate_pairs.py

Generates duplicate_pairs.json from experiments/results/duplicate_pairs.csv.
Captures metadata (source_file, source_commit, criteria, count) and
a bidirectional lookup map for fast O(1) query matching.
"""

from collections import defaultdict
import csv
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
CSV_PATH = REPO_ROOT / "experiments" / "results" / "duplicate_pairs.csv"
OUTPUT_JSON_PATH = REPO_ROOT / "duplicate_pairs.json"

SOURCE_COMMIT = "885e59a37ccb85139b43d54894d845177a895b34"
CRITERIA = [
    "generic_base_name",
    "first_block_inchikey",
    "desalted_canonical_smiles",
]


def generate_duplicate_pairs():
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Source CSV not found: {CSV_PATH}")

    pairs = []
    lookup = defaultdict(set)
    distinct_cids = set()

    with open(CSV_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cid1 = row["chembl_id_1"].strip()
            name1 = row["name_1"].strip()
            cid2 = row["chembl_id_2"].strip()
            name2 = row["name_2"].strip()
            raw_reasons = row["reasons"].strip()
            reasons = [r.strip() for r in raw_reasons.split(";") if r.strip()]

            pairs.append({
                "chembl_id_1": cid1,
                "name_1": name1,
                "chembl_id_2": cid2,
                "name_2": name2,
                "reasons": reasons,
            })

            lookup[cid1].add(cid2)
            lookup[cid2].add(cid1)
            distinct_cids.add(cid1)
            distinct_cids.add(cid2)

    count = len(pairs)
    sorted_lookup = {cid: sorted(list(dups)) for cid, dups in sorted(lookup.items())}

    data = {
        "metadata": {
            "source_file": "experiments/results/duplicate_pairs.csv",
            "source_commit": SOURCE_COMMIT,
            "criteria": CRITERIA,
            "count": count,
            "num_distinct_drugs": len(distinct_cids),
        },
        "pairs": pairs,
        "lookup": sorted_lookup,
    }

    with open(OUTPUT_JSON_PATH, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(f"Generated {OUTPUT_JSON_PATH}: {count} pairs across {len(distinct_cids)} distinct drugs.")
    return data


if __name__ == "__main__":
    generate_duplicate_pairs()
