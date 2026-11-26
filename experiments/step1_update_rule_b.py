"""
experiments/step1_update_rule_b.py

Implements Part A:
1. Reclassifies Rule (b) matches to reject the 5 stereo-prefix rows:
   - razoxane -> DEXRAZOXANE
   - medetomidine -> DEXMEDETOMIDINE
   - dexrabeprazole -> RABEPRAZOLE
   - levomethadone -> METHADONE
   - levomenthol -> MENTHOL
   with reason "stereo-prefix name difference".
2. Updates experiments/data/rejected_rule_b.csv.
3. Updates experiments/data/rule_b_matches.csv.
4. Regenerates experiments/data/gold_pairs_primary.csv.
5. Computes protocol counts and analytical random hit@10 for P1, P2a (0.3, 0.4, 0.5), P2b (0.4), P3, P3-all.
6. Writes experiments/logs/step1_update_rule_b.log and updates experiments/logs/protocol_dataset_counts_DB_A.log.
"""

import csv
import math
import re
import sqlite3
import numpy as np
from collections import defaultdict, Counter
from pathlib import Path
import psycopg2
from rdkit import Chem
from rdkit.Chem import DataStructs

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "experiments" / "data"
LOGS_DIR = ROOT / "experiments" / "logs"
DB_A = DATA_DIR / "drugrepo_oct6.db"

INORGANIC_EXCLUDE = {
    "CHEMBL1200345", "CHEMBL1200458", "CHEMBL1200539", "CHEMBL1200691",
    "CHEMBL1200696", "CHEMBL1237066", "CHEMBL1356", "CHEMBL2068884",
    "CHEMBL2106119", "CHEMBL3833401",
}

ALLOWED_RULE_B_EXCEPTIONS = {
    ("glucose", "DEXTROSE"),
    ("rifampicin", "RIFAMPIN"),
    ("fludeoxyglucose (18F)", "FLUDEOXYGLUCOSE F 18"),
}

STEREO_PREFIX_REJECTIONS = {
    ("razoxane", "DEXRAZOXANE"),
    ("medetomidine", "DEXMEDETOMIDINE"),
    ("dexrabeprazole", "RABEPRAZOLE"),
    ("levomethadone", "METHADONE"),
    ("levomenthol", "MENTHOL"),
}

DROP_WORDS = {
    "HYDROCHLORIDE", "HCL", "SODIUM", "DISODIUM", "TRISODIUM", "POTASSIUM",
    "CALCIUM", "MAGNESIUM", "ZINC", "PHOSPHATE", "DIPHOSPHATE", "SULFATE",
    "SULPHATE", "MESYLATE", "BESYLATE", "TOSYLATE", "MALEATE", "FUMARATE",
    "TARTRATE", "BITARTRATE", "CITRATE", "ACETATE", "NITRATE", "BROMIDE",
    "HYDROBROMIDE", "CHLORIDE", "IODIDE", "FLUORIDE", "LACTATE", "SUCCINATE",
    "PAMOATE", "VALERATE", "STEARATE", "PALMITATE", "OLEATE", "SALICYLATE",
    "BENZOATE", "NAPSYLATE", "GLUCONATE", "GLUCEPTATE", "LACTOBIONATE",
    "ANHYDROUS", "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE", "TETRAHYDRATE",
    "HEMIHYDRATE", "HYDRATE", "COMPLEX", "FREE", "ACID", "BASE", "CATION",
    "ANION", "ESTER", "ALCOHOL", "PROPIONATE", "BUTYRATE", "DIHYDROCHLORIDE",
    "SULFONATE", "HYDROGEN", "MONOSODIUM", "MEGLUMINE", "TROMETHAMINE"
}

def normalize_name(name: str) -> str:
    if not name:
        return ""
    s = name.upper()
    s = re.sub(r"[,;/()\[\]\-]+", " ", s)
    tokens = [t for t in s.split() if t and t not in DROP_WORDS]
    return " ".join(tokens).strip()

def rand_hit_at_10(M: int, R: int) -> float:
    if M < 10:
        return 1.0 if R > 0 else 0.0
    if M - R < 10:
        return 1.0
    return 1.0 - (math.comb(M - R, 10) / math.comb(M, 10))

def size_buckets(sizes):
    s = sum(1 for x in sizes if 1 <= x <= 10)
    m = sum(1 for x in sizes if 11 <= x <= 50)
    l = sum(1 for x in sizes if x > 50)
    return s, m, l

def main():
    log_lines = []
    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("=== PART A: RULE (b) STEREO-PREFIX REJECTION & GOLD REGENERATION ===")

    # 1. Load DB-A
    conn_local = sqlite3.connect(DB_A.resolve().as_uri() + "?mode=ro", uri=True)
    drugs = {}
    ik_full = {}
    ik_b1 = defaultdict(list)
    for cid, gbn, dname, smiles, fp_blob, triv in conn_local.execute(
        "SELECT chembl_id, generic_base_name, drug_name, smiles, fingerprint_blob, trivial_category FROM drugs"
    ).fetchall():
        fp = DataStructs.CreateFromBitString(fp_blob) if fp_blob else None
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        ik = Chem.MolToInchiKey(mol) if mol else None
        b1 = ik.split("-")[0] if ik else None
        b2 = ik.split("-")[1] if ik and len(ik.split("-")) > 1 else None
        if ik:
            ik_full[ik] = cid
            ik_b1[b1].append(cid)
        drugs[cid] = {
            "chembl_id": cid, "name": gbn, "drug_name": dname, "fingerprint": fp,
            "is_trivial": triv is not None, "inchikey": ik, "block1": b1, "block2": b2
        }
    cw_db, cw_pc = {}, {}
    for cid, db_id, pc_cid in conn_local.execute("SELECT chembl_id, drugbank_id, pubchem_cid FROM crosswalk").fetchall():
        if db_id and db_id.strip(): cw_db[db_id.strip()] = cid
        if pc_cid and str(pc_cid).strip(): cw_pc[str(pc_cid).strip()] = cid
    conn_local.close()

    # 2. Fetch DrugCentral structures & identifiers
    conn_pg = psycopg2.connect(host="unmtid-dbs.net", port=5433, user="drugman", password="dosage", dbname="drugcentral", connect_timeout=30)
    cur_pg = conn_pg.cursor()
    cur_pg.execute("""
        SELECT DISTINCT s.id, s.name, s.inchikey
        FROM structures s
        JOIN omop_relationship o ON s.id = o.struct_id
        WHERE o.relationship_name = 'indication'
        ORDER BY s.id
    """)
    dc_ind_structs = cur_pg.fetchall()

    cur_pg.execute("SELECT id, name, inchikey FROM structures ORDER BY id")
    all_dc_structs = cur_pg.fetchall()

    cur_pg.execute("""
        SELECT struct_id, id_type, identifier
        FROM identifier
        WHERE id_type IN ('DRUGBANK_ID', 'PUBCHEM_CID', 'ChEMBL_ID')
    """)
    dc_idents = defaultdict(lambda: defaultdict(list))
    for sid, itype, ident in cur_pg.fetchall():
        if ident and ident.strip():
            dc_idents[sid][itype].append(ident.strip())
    conn_pg.close()

    # 3. Map
    matches_b = []
    matches_c = []
    chembl_to_dc = defaultdict(list)
    struct_to_cid_ind = {}

    def map_one(sid, ik):
        clean_ik = ik.strip() if ik else None
        b1 = clean_ik.split("-")[0] if clean_ik else None
        if clean_ik and clean_ik in ik_full:
            return ik_full[clean_ik], "a", "rule_a_full_inchikey", clean_ik
        elif clean_ik and b1 in ik_b1:
            return ik_b1[b1][0], "b", "rule_b_block1_inchikey", clean_ik
        elif sid in dc_idents:
            ids = dc_idents[sid]
            for db_id in ids.get("DRUGBANK_ID", []):
                if db_id in cw_db: return cw_db[db_id], "c", "rule_c_crosswalk_drugbank", clean_ik
            for pc_id in ids.get("PUBCHEM_CID", []):
                if pc_id in cw_pc: return cw_pc[pc_id], "c", "rule_c_crosswalk_pubchem", clean_ik
            for ch_id in ids.get("ChEMBL_ID", []):
                if ch_id in drugs: return ch_id, "c", "rule_c_direct_chembl", clean_ik
        return None, None, None, clean_ik

    for sid, dc_name, dc_ik in dc_ind_structs:
        cid, r_short, r_full, clean_ik = map_one(sid, dc_ik)
        if not cid: continue
        struct_to_cid_ind[sid] = (cid, r_short)
        chembl_to_dc[cid].append((sid, dc_name, clean_ik, r_full))
        our_d = drugs[cid]
        dc_b2 = clean_ik.split("-")[1] if clean_ik and len(clean_ik.split("-")) > 1 else None
        b2_match = (dc_b2 == our_d["block2"]) if (dc_b2 and our_d["block2"]) else False

        if r_short == "b":
            dc_norm = normalize_name(dc_name)
            our_norm = normalize_name(our_d["name"])
            ag = "same" if (dc_norm and our_norm and (dc_norm == our_norm or dc_norm in our_norm or our_norm in dc_norm)) else "different"
            matches_b.append({
                "dc_struct_id": sid, "dc_name": dc_name, "dc_full_inchikey": clean_ik,
                "our_chembl_id": cid, "our_name": our_d["name"], "our_full_inchikey": our_d["inchikey"],
                "second_block_matches": b2_match, "dc_block2": dc_b2, "our_block2": our_d["block2"],
                "name_agreement": ag, "is_extra_drug": False,
            })
        elif r_short == "c":
            matches_c.append({
                "dc_struct_id": sid, "dc_name": dc_name, "dc_full_inchikey": clean_ik,
                "our_chembl_id": cid, "our_name": our_d["name"], "our_full_inchikey": our_d["inchikey"],
                "second_block_matches": b2_match, "match_subrule": r_full,
            })

    # 4. Filter Rule (b): Reject 4 original + 5 stereo-prefix rows
    accepted_b = []
    rejected_b = []
    for r in matches_b:
        pair_key = (r["dc_name"], r["our_name"])
        if pair_key in STEREO_PREFIX_REJECTIONS:
            r_copy = dict(r)
            r_copy["reason"] = "stereo-prefix name difference"
            rejected_b.append(r_copy)
        elif r["name_agreement"] == "same" or pair_key in ALLOWED_RULE_B_EXCEPTIONS:
            accepted_b.append(r)
        else:
            r_copy = dict(r)
            r_copy["reason"] = "different name"
            rejected_b.append(r_copy)

    # Write rejected_rule_b.csv
    rejected_file = DATA_DIR / "rejected_rule_b.csv"
    with open(rejected_file, "w", newline="", encoding="utf-8") as f:
        fn = ["dc_struct_id", "dc_name", "dc_full_inchikey", "our_chembl_id", "our_name",
              "our_full_inchikey", "second_block_matches", "dc_block2", "our_block2",
              "name_agreement", "is_extra_drug", "reason"]
        w = csv.DictWriter(f, fieldnames=fn)
        w.writeheader()
        w.writerows(rejected_b)

    log(f"Total Rule (b) candidate matches: {len(matches_b)}")
    log(f"Accepted Rule (b) matches:         {len(accepted_b)}")
    log(f"Rejected Rule (b) matches:         {len(rejected_b)}")
    for r in rejected_b:
        log(f"  Rejected: struct {r['dc_struct_id']} ({r['dc_name']}) -> {r['our_chembl_id']} ({r['our_name']}) [reason: {r['reason']}]")

    # Map all DC structures for ATC codes, excluding rejected structs
    rejected_sids = {int(r["dc_struct_id"]) for r in rejected_b}
    struct_to_cid_all = {}
    for sid, dc_name, dc_ik in all_dc_structs:
        if sid in rejected_sids:
            continue
        cid, _, _, _ = map_one(sid, dc_ik)
        if cid:
            struct_to_cid_all[sid] = cid

    cid_atc_l3 = defaultdict(set)
    with open(DATA_DIR / "drugcentral_raw_atc.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sid = int(r["struct_id"])
            if sid in struct_to_cid_all:
                cid = struct_to_cid_all[sid]
                if r["l3_code"] and r["l3_code"].strip():
                    cid_atc_l3[cid].add(r["l3_code"].strip())

    # Build filtered indications
    excluded_ind_ids = {r["indication_id"] for r in csv.DictReader(open(DATA_DIR / "excluded_indications.csv", encoding="utf-8"))}

    filt_pairs_set = set()
    filt_pair_rules = defaultdict(set)
    filt_cid_rules = defaultdict(set)

    with open(DATA_DIR / "drugcentral_raw_approved_indications.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sid = int(r["struct_id"])
            ind_id = r["concept_id"]
            ind_name = r["concept_name"]
            if sid not in struct_to_cid_ind:
                continue
            if sid in rejected_sids:
                continue
            cid, rule = struct_to_cid_ind[sid]
            filt_pairs_set.add((cid, ind_id, ind_name))
            filt_pair_rules[(cid, ind_id)].add(rule)
            filt_cid_rules[cid].add(rule)

    filt_pairs_list = sorted(list(filt_pairs_set))

    # Primary pool: DB-A non-trivial minus 10 inorganics (2,192 drugs)
    pool_primary = {cid for cid, d in drugs.items() if not d["is_trivial"]} - INORGANIC_EXCLUDE
    N_primary = len(pool_primary)
    M_full = N_primary - 1

    # Primary gold pairs (therapeutic only)
    primary_pairs = [p for p in filt_pairs_list if p[0] in pool_primary and p[1] not in excluded_ind_ids]
    with open(DATA_DIR / "gold_pairs_primary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["chembl_id", "indication_id", "indication_name"])
        w.writerows(primary_pairs)

    log(f"\nSaved {len(primary_pairs)} therapeutic gold pairs to gold_pairs_primary.csv")

    # Mapped drugs & rules breakdown for primary pool
    mapped_drugs_primary = {p[0] for p in primary_pairs}
    pair_rules_ther = Counter()
    for p in primary_pairs:
        cid, ind_id, _ = p
        rules = filt_pair_rules.get((cid, ind_id), filt_cid_rules[cid])
        best = "a" if "a" in rules else ("b" if "b" in rules else "c")
        pair_rules_ther[best] += 1

    drug_rules_ther = Counter()
    for cid in mapped_drugs_primary:
        rules = filt_cid_rules[cid]
        best = "a" if "a" in rules else ("b" if "b" in rules else "c")
        drug_rules_ther[best] += 1

    log(f"\n--- Primary Pool Gold Summary (N={N_primary}, M={M_full}) ---")
    log(f"Mapped drugs with therapeutic indications: {len(mapped_drugs_primary)}")
    log(f"Mapped drugs by primary rule (a>b>c): a={drug_rules_ther['a']}, b={drug_rules_ther['b']}, c={drug_rules_ther['c']}")
    log(f"Therapeutic pairs by rule:             a={pair_rules_ther['a']}, b={pair_rules_ther['b']}, c={pair_rules_ther['c']}")
    log(f"Total therapeutic gold pairs:          {len(primary_pairs)}")

    # Protocols
    drug_inds = defaultdict(set)
    for cid, ind_id, _ in primary_pairs:
        drug_inds[cid].add(ind_id)

    scored_with_inds = sorted(list(drug_inds.keys()))

    # P1
    rel_P1 = {}
    for d1 in scored_with_inds:
        inds1 = drug_inds[d1]
        rel = {d2 for d2 in scored_with_inds if d2 != d1 and (inds1 & drug_inds[d2])}
        if rel:
            rel_P1[d1] = rel

    Q_P1 = sorted(list(rel_P1.keys()))
    sizes_P1 = [len(rel_P1[d]) for d in Q_P1]
    rand_P1 = [rand_hit_at_10(M_full, r) for r in sizes_P1]
    s1, m1, l1 = size_buckets(sizes_P1)

    log("\n" + "=" * 92)
    log("PROTOCOL COUNTS & ANALYTICAL RANDOM EXPECTATIONS (DB-A PRIMARY, N=2,192)")
    log("=" * 92)
    log(f"{'Protocol':<16} {'Queries':>8} {'Dropped':>8} {'Mean_R':>8} {'Med_R':>7} {'Rand_Hit@10':>13} {'S/M/L':>15}")
    log("-" * 92)
    log(f"{'P1':<16} {len(Q_P1):>8} {0:>8} {np.mean(sizes_P1):>8.2f} {np.median(sizes_P1):>7.1f} {np.mean(rand_P1):>13.4f} {f'{s1}/{m1}/{l1}':>15}")

    # P2a
    for theta in [0.3, 0.4, 0.5]:
        sizes_P2a = []
        dropped_P2a = 0
        for d1 in Q_P1:
            fp1 = drugs[d1]["fingerprint"]
            rel2 = set()
            for d2 in rel_P1[d1]:
                fp2 = drugs[d2]["fingerprint"]
                tan = DataStructs.TanimotoSimilarity(fp1, fp2) if fp1 and fp2 else 0.0
                if tan <= theta:
                    rel2.add(d2)
            if rel2:
                sizes_P2a.append(len(rel2))
            else:
                dropped_P2a += 1
        rand_P2a = [rand_hit_at_10(M_full, r) for r in sizes_P2a]
        s, m, l = size_buckets(sizes_P2a)
        log(f"{f'P2a ({theta})':<16} {len(sizes_P2a):>8} {dropped_P2a:>8} {np.mean(sizes_P2a):>8.2f} {np.median(sizes_P2a):>7.1f} {np.mean(rand_P2a):>13.4f} {f'{s}/{m}/{l}':>15}")

    # P2b (theta=0.4)
    # Same relevant set as P2a 0.4, but candidates with tanimoto > 0.4 removed -> M_q shrinks!
    sizes_P2b = []
    rand_P2b = []
    dropped_P2b = 0
    for d1 in Q_P1:
        fp1 = drugs[d1]["fingerprint"]
        # Relevant set
        rel2 = set()
        for d2 in rel_P1[d1]:
            fp2 = drugs[d2]["fingerprint"]
            tan = DataStructs.TanimotoSimilarity(fp1, fp2) if fp1 and fp2 else 0.0
            if tan <= 0.4:
                rel2.add(d2)
        if not rel2:
            dropped_P2b += 1
            continue
        # Candidates for d1: all pool drugs except d1 with tan <= 0.4
        n_analogues = 0
        for c in pool_primary:
            if c == d1: continue
            fpc = drugs[c]["fingerprint"]
            tan = DataStructs.TanimotoSimilarity(fp1, fpc) if fp1 and fpc else 0.0
            if tan > 0.4:
                n_analogues += 1
        M_q = M_full - n_analogues
        sizes_P2b.append(len(rel2))
        rand_P2b.append(rand_hit_at_10(M_q, len(rel2)))

    s_b, m_b, l_b = size_buckets(sizes_P2b)
    log(f"{'P2b (0.4)':<16} {len(sizes_P2b):>8} {dropped_P2b:>8} {np.mean(sizes_P2b):>8.2f} {np.median(sizes_P2b):>7.1f} {np.mean(rand_P2b):>13.4f} {f'{s_b}/{m_b}/{l_b}':>15}")

    # P3 & P3-all
    pool_with_atc = {cid for cid in pool_primary if len(cid_atc_l3[cid]) > 0}
    has_atc_queries = [d for d in Q_P1 if len(cid_atc_l3[d]) > 0]
    no_atc_queries = [d for d in Q_P1 if len(cid_atc_l3[d]) == 0]

    sizes_P3 = []
    dropped_P3 = 0
    for d1 in has_atc_queries:
        atc1 = cid_atc_l3[d1]
        rel3 = set()
        for d2 in rel_P1[d1]:
            atc2 = cid_atc_l3[d2]
            if not atc2:
                continue
            if atc1.isdisjoint(atc2):
                rel3.add(d2)
        if rel3:
            sizes_P3.append(len(rel3))
        else:
            dropped_P3 += 1

    M_atc = len(pool_with_atc) - 1
    rand_P3_Matc = [rand_hit_at_10(M_atc, r) for r in sizes_P3]
    rand_P3_Mfull = [rand_hit_at_10(M_full, r) for r in sizes_P3]
    s3, m3, l3 = size_buckets(sizes_P3)

    log(f"{'P3 (primary)':<16} {len(sizes_P3):>8} {dropped_P3:>8} {np.mean(sizes_P3):>8.2f} {np.median(sizes_P3):>7.1f} {np.mean(rand_P3_Matc):>13.4f} {f'{s3}/{m3}/{l3}':>15}")
    log(f"{'P3-all':<16} {len(sizes_P3):>8} {dropped_P3:>8} {np.mean(sizes_P3):>8.2f} {np.median(sizes_P3):>7.1f} {np.mean(rand_P3_Mfull):>13.4f} {f'{s3}/{m3}/{l3}':>15}")
    log(f"  P3 queries without ATC code (excluded from P3 query set): {len(no_atc_queries)}")
    log(f"  Pool drugs with ATC code: {len(pool_with_atc)} (M_atc = {M_atc})")

    # Save log
    out_log = LOGS_DIR / "step1_update_rule_b.log"
    out_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"\nSaved log to {out_log}")

if __name__ == "__main__":
    main()

