"""
experiments/step1_gold_standard.py — DrugCentral Gold Standard Builder

Step 1 of the MoleculAR evaluation harness:
1. Extracts approved indications from DrugCentral (omop_relationship WHERE relationship_name = 'indication').
2. Maps DrugCentral drugs to chembl_id in priority order:
   (a) Full InChIKey match
   (b) First block (14-char) InChIKey match
   (c) External IDs via crosswalk table (DrugBank, PubChem, direct ChEMBL)
3. Saves:
   - experiments/data/gold_source.txt
   - experiments/data/gold_pairs.csv (chembl_id, indication_id, indication_name)
   - experiments/data/unmapped.csv
   - experiments/logs/step1_gold_standard.log
"""

import argparse
import csv
import hashlib
import random
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from collections import Counter

import psycopg2
from rdkit import Chem


# DrugCentral connection details from drugcentral.org/download
DC_HOST = "unmtid-dbs.net"
DC_PORT = 5433
DC_USER = "drugman"
DC_PASS = "dosage"
DC_DB   = "drugcentral"
DC_URL  = "https://drugcentral.org/download (PostgreSQL instance drugcentral:unmtid-dbs.net:5433)"


def get_sha256(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def run_step1(db_path: str, seed: int = 42):
    random.seed(seed)
    data_dir = Path(db_path).parent
    logs_dir = data_dir.parent / "logs"
    data_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_file = logs_dir / "step1_gold_standard.log"
    log_fp = open(log_file, "w", encoding="utf-8")

    def log(msg: str):
        print(msg)
        log_fp.write(msg + "\n")
        log_fp.flush()

    log(f"=== STEP 1: GOLD STANDARD EXTRACTION AND MAPPING ===")
    log(f"Timestamp: {datetime.now(timezone.utc).isoformat()}")
    log(f"Random seed: {seed}")
    log(f"Local database snapshot: {db_path}")

    # -------------------------------------------------------------
    # 1. Load Local Snapshot Database
    # -------------------------------------------------------------
    log("\n--- 1. Loading local database snapshot ---")
    conn_local = sqlite3.connect(Path(db_path).resolve().as_uri() + "?mode=ro", uri=True)
    cur_local = conn_local.cursor()

    cur_local.execute("SELECT chembl_id, smiles, generic_base_name, trivial_category FROM drugs")
    local_drugs = {}
    local_ik_full = {}
    local_ik_block1 = {}

    for cid, smiles, gbn, triv in cur_local.fetchall():
        is_trivial = triv is not None
        ik = None
        b1 = None
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                ik = Chem.MolToInchiKey(mol)
                if ik:
                    b1 = ik.split("-")[0]
                    local_ik_full[ik] = cid
                    local_ik_block1.setdefault(b1, []).append(cid)

        local_drugs[cid] = {
            "chembl_id": cid,
            "generic_base_name": gbn,
            "is_trivial": is_trivial,
            "trivial_category": triv,
            "inchikey": ik,
            "block1": b1,
        }

    # Load crosswalk mappings
    cur_local.execute("SELECT chembl_id, drugbank_id, pubchem_cid FROM crosswalk")
    cw_drugbank = {}
    cw_pubchem = {}
    for cid, db_id, pc_cid in cur_local.fetchall():
        if db_id and db_id.strip():
            cw_drugbank[db_id.strip()] = cid
        if pc_cid and str(pc_cid).strip():
            cw_pubchem[str(pc_cid).strip()] = cid

    conn_local.close()

    total_local = len(local_drugs)
    scored_cids = {cid for cid, d in local_drugs.items() if not d["is_trivial"]}
    trivial_cids = {cid for cid, d in local_drugs.items() if d["is_trivial"]}
    log(f"Local total drugs: {total_local}")
    log(f"Local non-trivial (scored pool): {len(scored_cids)}")
    log(f"Local trivial: {len(trivial_cids)}")
    log(f"Local with full InChIKey: {len(local_ik_full)}")
    log(f"Local unique Block1: {len(local_ik_block1)}")
    log(f"Crosswalk DrugBank mappings: {len(cw_drugbank)}")
    log(f"Crosswalk PubChem mappings: {len(cw_pubchem)}")

    # -------------------------------------------------------------
    # 2. Connect to DrugCentral PostgreSQL
    # -------------------------------------------------------------
    log("\n--- 2. Connecting to DrugCentral PostgreSQL ---")
    log(f"Host: {DC_HOST}:{DC_PORT}, Database: {DC_DB}, User: {DC_USER}")
    try:
        conn_pg = psycopg2.connect(
            host=DC_HOST,
            port=DC_PORT,
            user=DC_USER,
            password=DC_PASS,
            dbname=DC_DB,
            connect_timeout=30,
        )
    except Exception as e:
        log(f"ERROR: Could not connect to DrugCentral: {e}")
        log_fp.close()
        sys.exit(1)

    cur_pg = conn_pg.cursor()

    # Save raw data export for reproducibility and SHA-256 tracking
    raw_ind_file = data_dir / "drugcentral_raw_approved_indications.csv"
    log(f"Dumping approved indications from omop_relationship (relationship_name = 'indication') to {raw_ind_file.name} ...")
    
    cur_pg.execute("""
        SELECT struct_id, concept_id, concept_name, umls_cui, snomed_conceptid, relationship_name
        FROM omop_relationship
        WHERE relationship_name = 'indication'
        ORDER BY struct_id, concept_id
    """)
    raw_ind_rows = cur_pg.fetchall()
    
    with open(raw_ind_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["struct_id", "concept_id", "concept_name", "umls_cui", "snomed_conceptid", "relationship_name"])
        writer.writerows(raw_ind_rows)

    raw_ind_sha256 = get_sha256(raw_ind_file)
    log(f"Raw approved indications count: {len(raw_ind_rows)}")
    log(f"Raw indications CSV SHA-256: {raw_ind_sha256}")

    # Fetch all DrugCentral structures
    cur_pg.execute("SELECT id, name, inchikey, smiles FROM structures ORDER BY id")
    all_dc_structures = cur_pg.fetchall()
    log(f"Total DrugCentral structures: {len(all_dc_structures)}")

    # Fetch external identifiers for mapping rule (c)
    cur_pg.execute("""
        SELECT struct_id, id_type, identifier
        FROM identifier
        WHERE id_type IN ('DRUGBANK_ID', 'PUBCHEM_CID', 'ChEMBL_ID')
        ORDER BY struct_id, id_type
    """)
    dc_identifiers = {}
    for sid, id_type, ident in cur_pg.fetchall():
        if ident and ident.strip():
            dc_identifiers.setdefault(sid, {}).setdefault(id_type, []).append(ident.strip())

    conn_pg.close()

    # Write gold_source.txt
    source_txt = data_dir / "gold_source.txt"
    with open(source_txt, "w", encoding="utf-8") as f:
        f.write("DRUGCENTRAL GOLD STANDARD SOURCE METADATA\n")
        f.write("=========================================\n")
        f.write("Website Download Portal: https://drugcentral.org/download\n")
        f.write("Published Release Name: DrugCentral 2027 Update (Human & Veterinary Drugs & Uses)\n")
        f.write("Published Database Dump Reference: Drugcentral_2026-09-25.pgdump (or drugcentral.dump.09022026.sql.gz)\n")
        f.write("Published Dump URL: https://unmtid-dbs.net/download/DrugCentral/2027/db_dump/Drugcentral_2026-09-25.pgdump\n")
        f.write("Published Dump Date: 2026-09-25\n")
        f.write("Database Engine: PostgreSQL 16.15\n")
        f.write(f"Queried Endpoint: drugcentral:{DC_HOST}:{DC_PORT} (user: \"{DC_USER}\", dbname: \"{DC_DB}\")\n")
        f.write("Endpoint Provenance: Explicitly published on https://drugcentral.org/download\n")
        f.write("Live Instance Confirmation: Official live public instance corresponding to DrugCentral 2027 release documented on drugcentral.org/download (bitwise comparison against 1.3 GB .pgdump not performed).\n")
        f.write(f"Download/Access Date: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}\n")
        f.write("Filter Column: omop_relationship.relationship_name\n")
        f.write("Filter Condition: relationship_name = 'indication' (approved indications only; excludes 'off-label use', 'contraindication', etc.)\n")
        f.write(f"Extracted Table: omop_relationship ({len(raw_ind_rows)} approved indication rows across {len({r[0] for r in raw_ind_rows})} distinct structures)\n")
        f.write(f"Raw Extracted CSV: {raw_ind_file.name}\n")
        f.write(f"Raw Extracted CSV SHA-256: {raw_ind_sha256}\n")

    log(f"Saved source metadata to {source_txt}")

    # -------------------------------------------------------------
    # 3. DrugCentral Drugs Mapping
    # -------------------------------------------------------------
    log("\n--- 3. Mapping DrugCentral drugs to chembl_id ---")

    # Distinct struct_ids that actually have approved indications
    dc_struct_ids_with_ind = {r[0] for r in raw_ind_rows}
    log(f"DrugCentral structures total: {len(all_dc_structures)}")
    log(f"DrugCentral structures with approved indications: {len(dc_struct_ids_with_ind)}")

    mapped_drugs = {}   # struct_id -> {'chembl_id', 'rule', 'name', 'inchikey'}
    unmapped_drugs = [] # list of dicts

    rule_counter = Counter()

    for sid, name, ik, smiles in all_dc_structures:
        has_ind = sid in dc_struct_ids_with_ind
        clean_ik = ik.strip() if ik else None
        matched_cid = None
        matched_rule = None

        # Rule (a): InChIKey full match
        if clean_ik and clean_ik in local_ik_full:
            matched_cid = local_ik_full[clean_ik]
            matched_rule = "rule_a_full_inchikey"

        # Rule (b): InChIKey first block (14 chars) match
        elif clean_ik:
            b1 = clean_ik.split("-")[0]
            if b1 in local_ik_block1:
                # If multiple local drugs share Block1, pick the first
                matched_cid = local_ik_block1[b1][0]
                matched_rule = "rule_b_block1_inchikey"

        # Rule (c): Crosswalk table external IDs
        if not matched_cid and sid in dc_identifiers:
            ids = dc_identifiers[sid]
            # c1: DrugBank
            for db_id in ids.get("DRUGBANK_ID", []):
                if db_id in cw_drugbank:
                    matched_cid = cw_drugbank[db_id]
                    matched_rule = "rule_c_crosswalk_drugbank"
                    break
            # c2: PubChem
            if not matched_cid:
                for pc_id in ids.get("PUBCHEM_CID", []):
                    if pc_id in cw_pubchem:
                        matched_cid = cw_pubchem[pc_id]
                        matched_rule = "rule_c_crosswalk_pubchem"
                        break
            # c3: Direct ChEMBL ID
            if not matched_cid:
                for ch_id in ids.get("ChEMBL_ID", []):
                    if ch_id in local_drugs:
                        matched_cid = ch_id
                        matched_rule = "rule_c_direct_chembl"
                        break

        if matched_cid:
            mapped_drugs[sid] = {
                "struct_id": sid,
                "chembl_id": matched_cid,
                "rule": matched_rule,
                "name": name,
                "inchikey": clean_ik,
                "has_indications": has_ind,
            }
            if has_ind:
                rule_counter[matched_rule] += 1
        else:
            unmapped_drugs.append({
                "struct_id": sid,
                "name": name,
                "inchikey": clean_ik,
                "has_indications": has_ind,
            })

    # Stats for DC drugs WITH indications
    mapped_with_ind = {sid: d for sid, d in mapped_drugs.items() if d["has_indications"]}
    unmapped_with_ind = [d for d in unmapped_drugs if d["has_indications"]]

    log(f"\nMapping Results (for DrugCentral drugs with approved indications):")
    log(f"  DrugCentral drugs with indications total: {len(dc_struct_ids_with_ind)}")
    log(f"  Mapped: {len(mapped_with_ind)}")
    log(f"  Unmapped: {len(unmapped_with_ind)}")
    log(f"  Mapping Rule Breakdown:")
    for rule, count in sorted(rule_counter.items()):
        log(f"    {rule}: {count}")

    # Check mapped drugs in our scored set
    mapped_scored = {d["chembl_id"] for d in mapped_with_ind.values() if d["chembl_id"] in scored_cids}
    mapped_trivial = {d["chembl_id"] for d in mapped_with_ind.values() if d["chembl_id"] in trivial_cids}
    log(f"\nLocal Drug Coverage (drugs with approved indications):")
    log(f"  Mapped drugs in our scored set (non-trivial): {len(mapped_scored)} unique chembl_ids")
    log(f"  Mapped drugs in trivial set: {len(mapped_trivial)} unique chembl_ids")

    # Also stats across ALL DrugCentral structures
    log(f"\nMapping Results (across ALL {len(all_dc_structures)} DrugCentral structures):")
    log(f"  Total mapped structures: {len(mapped_drugs)}")
    log(f"  Total unmapped structures: {len(unmapped_drugs)}")

    # -------------------------------------------------------------
    # 4. Save unmapped.csv
    # -------------------------------------------------------------
    unmapped_file = data_dir / "unmapped.csv"
    with open(unmapped_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["struct_id", "drug_name", "inchikey", "has_approved_indications"])
        for d in unmapped_drugs:
            writer.writerow([d["struct_id"], d["name"], d["inchikey"], d["has_indications"]])

    log(f"\nSaved unmapped structures to {unmapped_file.name} (SHA-256: {get_sha256(unmapped_file)})")

    # -------------------------------------------------------------
    # 5. Build and Save gold_pairs.csv
    # -------------------------------------------------------------
    # Schema requested: (chembl_id, indication_id, indication_name)
    log("\n--- 4. Building gold_pairs.csv ---")
    
    # We build pairs for mapped drugs. Note: per R8, we keep all mapped drugs or scored?
    # The prompt states:
    # "Save experiments/data/gold_pairs.csv (chembl_id, indication_id, indication_name) and experiments/data/unmapped.csv."
    # We save distinct (chembl_id, indication_id, indication_name) pairs.
    gold_pairs_set = set()
    scored_pairs_set = set()

    for r in raw_ind_rows:
        sid = r[0]
        concept_id = r[1]
        concept_name = r[2]
        if sid in mapped_drugs:
            cid = mapped_drugs[sid]["chembl_id"]
            pair = (cid, concept_id, concept_name)
            gold_pairs_set.add(pair)
            if cid in scored_cids:
                scored_pairs_set.add(pair)

    gold_pairs_list = sorted(list(gold_pairs_set), key=lambda x: (x[0], x[1]))
    
    gold_pairs_file = data_dir / "gold_pairs.csv"
    with open(gold_pairs_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["chembl_id", "indication_id", "indication_name"])
        for p in gold_pairs_list:
            writer.writerow(p)

    gold_pairs_sha256 = get_sha256(gold_pairs_file)
    log(f"Saved {len(gold_pairs_list)} distinct gold pairs to {gold_pairs_file.name}")
    log(f"gold_pairs.csv SHA-256: {gold_pairs_sha256}")

    # Indication statistics
    # 1. Total distinct indications across all mapped pairs
    distinct_indications_all = {p[1] for p in gold_pairs_set}
    counts_per_ind_all = Counter(p[1] for p in gold_pairs_set)
    inds_ge_3_all = sum(1 for c, cnt in counts_per_ind_all.items() if cnt >= 3)

    # 2. Total distinct indications in scored set
    distinct_indications_scored = {p[1] for p in scored_pairs_set}
    counts_per_ind_scored = Counter(p[1] for p in scored_pairs_set)
    inds_ge_3_scored = sum(1 for c, cnt in counts_per_ind_scored.items() if cnt >= 3)
    inds_ge_2_scored = sum(1 for c, cnt in counts_per_ind_scored.items() if cnt >= 2)

    log(f"\n--- Indication Statistics ---")
    log(f"In Scored Set (non-trivial drugs only):")
    log(f"  Distinct indications: {len(distinct_indications_scored)}")
    log(f"  Indications with >= 3 mapped drugs: {inds_ge_3_scored}")
    log(f"  Indications with >= 2 mapped drugs: {inds_ge_2_scored}")
    log(f"  Scored drugs with >= 1 indication: {len(mapped_scored)}")
    log(f"Across ALL Mapped Drugs (including trivial):")
    log(f"  Distinct indications: {len(distinct_indications_all)}")
    log(f"  Indications with >= 3 mapped drugs: {inds_ge_3_all}")

    # -------------------------------------------------------------
    # 6. Sample 10 Random Mapped Rows (seed 42)
    # -------------------------------------------------------------
    log(f"\n--- 10 Randomly Chosen Mapped Rows from gold_pairs.csv (seed={seed}) ---")
    rng = random.Random(seed)
    sample_indices = sorted(rng.sample(range(len(gold_pairs_list)), 10))
    sample_rows = [gold_pairs_list[i] for i in sample_indices]

    log(f"{'#':<3} {'chembl_id':<16} {'indication_id':<14} {'indication_name'}")
    log("-" * 75)
    for idx, (cid, ind_id, ind_name) in enumerate(sample_rows, 1):
        log(f"{idx:<3} {cid:<16} {ind_id:<14} {ind_name}")

    # Update manifest.json with gold standard details
    manifest_path = data_dir.parent / "results" / "manifest.json"
    if manifest_path.exists():
        import json
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest = json.load(f)
        manifest["gold_standard"] = {
            "source_url": DC_URL,
            "filter_column": "omop_relationship.relationship_name",
            "filter_value": "indication",
            "gold_pairs_csv": str(gold_pairs_file),
            "gold_pairs_sha256": gold_pairs_sha256,
            "raw_indications_csv": str(raw_ind_file),
            "raw_indications_sha256": raw_ind_sha256,
            "unmapped_csv": str(unmapped_file),
            "unmapped_sha256": get_sha256(unmapped_file),
            "drugs_with_indications_total": len(dc_struct_ids_with_ind),
            "drugs_with_indications_mapped": len(mapped_with_ind),
            "drugs_with_indications_unmapped": len(unmapped_with_ind),
            "mapped_in_scored_pool": len(mapped_scored),
            "distinct_indications_scored": len(distinct_indications_scored),
            "indications_ge_3_drugs_scored": inds_ge_3_scored,
        }
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)

    log(f"\nStep 1 complete. Manifest updated. Log written to {log_file}")
    log_fp.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default="experiments/data/drugrepo_snapshot.db", help="Path to database snapshot")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    run_step1(args.db, args.seed)

