import csv
import sqlite3
import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs
from collections import defaultdict
from pathlib import Path

data_dir = Path(__file__).resolve().parent / "data"
db_path = data_dir / "drugrepo_snapshot.db"

# 1. Load local drugs
conn_local = sqlite3.connect(db_path)
cur_local = conn_local.cursor()
cur_local.execute("SELECT chembl_id, generic_base_name, smiles, trivial_category, fingerprint_blob FROM drugs")

local_drugs = {}
local_fps = {}
local_ik_full = {}
local_ik_block1 = {}

for cid, gbn, smiles, triv, fp_blob in cur_local.fetchall():
    is_trivial = triv is not None
    fp = DataStructs.CreateFromBitString(fp_blob) if fp_blob else None
    local_drugs[cid] = {
        "chembl_id": cid,
        "name": gbn,
        "is_trivial": is_trivial,
    }
    local_fps[cid] = fp
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            ik = Chem.MolToInchiKey(mol)
            if ik:
                b1 = ik.split("-")[0]
                local_ik_full[ik] = cid
                local_ik_block1.setdefault(b1, []).append(cid)

cur_local.execute("SELECT chembl_id, drugbank_id, pubchem_cid FROM crosswalk")
cw_db = {}
cw_pc = {}
for cid, db_id, pc_cid in cur_local.fetchall():
    if db_id and db_id.strip(): cw_db[db_id.strip()] = cid
    if pc_cid and str(pc_cid).strip(): cw_pc[str(pc_cid).strip()] = cid
conn_local.close()

scored_cids = {cid for cid, d in local_drugs.items() if not d["is_trivial"]}

# 2. Build struct_to_cid mapping from rule files
struct_to_cid = {}

# From rule_b_matches.csv
with open(data_dir / "rule_b_matches.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for r in reader:
        struct_to_cid[int(r["dc_struct_id"])] = r["our_chembl_id"]

# From rule_c_matches.csv
with open(data_dir / "rule_c_matches.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for r in reader:
        struct_to_cid[int(r["dc_struct_id"])] = r["our_chembl_id"]

# From drugcentral_raw_approved_indications and structures for rule_a
# Let's read structures from postgres or query
# Actually, we can load all structures from drugcentral_raw_atc struct_ids
import psycopg2
conn_pg = psycopg2.connect(host="unmtid-dbs.net", port=5433, user="drugman", password="dosage", dbname="drugcentral", connect_timeout=30)
cur_pg = conn_pg.cursor()

cur_pg.execute("SELECT id, inchikey FROM structures WHERE inchikey IS NOT NULL")
for sid, ik in cur_pg.fetchall():
    clean_ik = ik.strip()
    if clean_ik in local_ik_full:
        struct_to_cid[sid] = local_ik_full[clean_ik]

# Load ATC codes from raw_atc_file
dc_atc_rows = []
with open(data_dir / "drugcentral_raw_atc.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for r in reader:
        dc_atc_rows.append((int(r["struct_id"]), r["l3_code"], r["atc_code"]))

conn_pg.close()

cid_atc_l3 = defaultdict(set)
for sid, l3, atc_code in dc_atc_rows:
    if sid in struct_to_cid:
        cid = struct_to_cid[sid]
        if l3 and l3.strip():
            cid_atc_l3[cid].add(l3.strip())

# 3. Load gold_pairs.csv (chembl_id, indication_id, indication_name)
# Exclude Diagnostic indications
excluded_ind_ids = set()
with open(data_dir / "excluded_indications.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for r in reader:
        excluded_ind_ids.add(r["indication_id"])

drug_indications = defaultdict(set)
with open(data_dir / "gold_pairs.csv", "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for r in reader:
        cid = r["chembl_id"]
        ind_id = r["indication_id"]
        if cid in scored_cids and ind_id not in excluded_ind_ids:
            drug_indications[cid].add(ind_id)

print(f"Scored drugs with therapeutic (non-diagnostic) indications: {len(drug_indications)}")
print(f"Scored drugs with at least one ATC L3 code: {len([cid for cid in scored_cids if len(cid_atc_l3[cid]) > 0])}")
print(f"Scored drugs with NO ATC code: {len([cid for cid in scored_cids if len(cid_atc_l3[cid]) == 0])}")

# 4. Build relevant sets for P1, P2 (0.3, 0.4, 0.5), P3
scored_list = sorted(list(drug_indications.keys()))
all_rel_P1 = {}

for d1 in scored_list:
    inds1 = drug_indications[d1]
    rel = set()
    for d2 in scored_list:
        if d1 == d2:
            continue
        if inds1 & drug_indications[d2]:
            rel.add(d2)
    if rel:
        all_rel_P1[d1] = rel

Q_P1 = sorted(list(all_rel_P1.keys()))
print(f"\n--- Protocol P1 ---")
print(f"Number of query drugs: {len(Q_P1)}")
sizes_P1 = [len(all_rel_P1[d]) for d in Q_P1]
print(f"P1 Mean rel size: {np.mean(sizes_P1):.2f}, Median rel size: {np.median(sizes_P1):.1f}")
print(f"P1 dropped queries (empty relevant set): 0")

# Protocol P2 with thresholds 0.3, 0.4, 0.5
p2_stats = {}
for theta in [0.3, 0.4, 0.5]:
    Q_P2 = []
    sizes_P2 = []
    dropped_P2 = 0
    for d1 in Q_P1:
        fp1 = local_fps[d1]
        rel_P2 = set()
        for d2 in all_rel_P1[d1]:
            fp2 = local_fps[d2]
            tan = DataStructs.TanimotoSimilarity(fp1, fp2) if fp1 and fp2 else 0.0
            if tan <= theta:
                rel_P2.add(d2)
        if len(rel_P2) > 0:
            Q_P2.append(d1)
            sizes_P2.append(len(rel_P2))
        else:
            dropped_P2 += 1
    p2_stats[theta] = (len(Q_P2), np.mean(sizes_P2), np.median(sizes_P2), dropped_P2)
    print(f"\n--- Protocol P2 (theta={theta}) ---")
    print(f"Number of query drugs: {len(Q_P2)}")
    print(f"Dropped queries (empty rel set): {dropped_P2}")
    print(f"Mean rel size: {np.mean(sizes_P2):.2f}, Median rel size: {np.median(sizes_P2):.1f}")

# Protocol P3:
# relevant drugs share an indication with query AND have a different ATC level-3 code from query
# shares NONE of the query's level-3 codes.
# Drugs with no ATC code are excluded from P3 and counted.
no_atc_in_Q = sum(1 for d in Q_P1 if len(cid_atc_l3[d]) == 0)
has_atc_Q = [d for d in Q_P1 if len(cid_atc_l3[d]) > 0]

Q_P3 = []
sizes_P3 = []
dropped_P3 = 0

for d1 in has_atc_Q:
    atc1 = cid_atc_l3[d1]
    rel_P3 = set()
    for d2 in all_rel_P1[d1]:
        atc2 = cid_atc_l3[d2]
        if len(atc2) == 0:
            continue
        if atc1.isdisjoint(atc2):
            rel_P3.add(d2)
    if len(rel_P3) > 0:
        Q_P3.append(d1)
        sizes_P3.append(len(rel_P3))
    else:
        dropped_P3 += 1

print(f"\n--- Protocol P3 ---")
print(f"Query drugs excluded from P3 due to NO ATC code: {no_atc_in_Q}")
print(f"Query drugs with ATC code evaluated: {len(has_atc_Q)}")
print(f"Number of query drugs in P3: {len(Q_P3)}")
print(f"Dropped queries (empty rel set among ATC-having): {dropped_P3}")
print(f"Total dropped from P1 baseline to P3: {no_atc_in_Q + dropped_P3}")
print(f"Mean rel size: {np.mean(sizes_P3):.2f}, Median rel size: {np.median(sizes_P3):.1f}")

# Save all these counts to a log file
log_path = data_dir.parent / "logs" / "protocol_dataset_counts.log"
with open(log_path, "w", encoding="utf-8") as f:
    f.write("DATASET COUNTS FOR PROTOCOLS P1, P2 (0.3, 0.4, 0.5), P3\n")
    f.write("========================================================\n")
    f.write(f"Scored drugs with therapeutic indications: {len(drug_indications)}\n")
    f.write(f"Scored drugs with at least one ATC L3 code: {len([cid for cid in scored_cids if len(cid_atc_l3[cid]) > 0])}\n")
    f.write(f"Scored drugs with NO ATC code: {len([cid for cid in scored_cids if len(cid_atc_l3[cid]) == 0])}\n\n")
    f.write(f"Protocol P1:\n  Query drugs: {len(Q_P1)}\n  Mean rel size: {np.mean(sizes_P1):.2f}\n  Median rel size: {np.median(sizes_P1):.1f}\n  Dropped queries: 0\n\n")
    for theta in [0.3, 0.4, 0.5]:
        q_len, m_size, med_size, drp = p2_stats[theta]
        f.write(f"Protocol P2 (theta={theta}):\n  Query drugs: {q_len}\n  Mean rel size: {m_size:.2f}\n  Median rel size: {med_size:.1f}\n  Dropped queries: {drp}\n\n")
    f.write(f"Protocol P3:\n  Queries excluded due to no ATC code: {no_atc_in_Q}\n  Query drugs: {len(Q_P3)}\n  Mean rel size: {np.mean(sizes_P3):.2f}\n  Median rel size: {np.median(sizes_P3):.1f}\n  Dropped queries (empty rel set): {dropped_P3}\n  Total dropped from P1: {no_atc_in_Q + dropped_P3}\n")
print(f"Saved: {log_path}")

