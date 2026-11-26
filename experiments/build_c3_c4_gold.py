import csv
import sqlite3
import psycopg2
from collections import defaultdict
from pathlib import Path

ROOT = Path(r"C:\Users\hg340\Documents\MoleculAR")
DATA_DIR = ROOT / "experiments" / "data"
DB_A = DATA_DIR / "drugrepo_oct6.db"

INORGANIC_EXCLUDE = {
    "CHEMBL1200345", "CHEMBL1200458", "CHEMBL1200539", "CHEMBL1200691",
    "CHEMBL1200696", "CHEMBL1237066", "CHEMBL1356", "CHEMBL2068884",
    "CHEMBL2106119", "CHEMBL3833401",
}

def main():
    conn_local = sqlite3.connect(DB_A.resolve().as_uri() + "?mode=ro", uri=True)
    ik_full = {}
    ik_b1 = defaultdict(list)
    drugs = {}
    for cid, smiles, triv in conn_local.execute("SELECT chembl_id, smiles, trivial_category FROM drugs").fetchall():
        drugs[cid] = {"triv": triv is not None}
        if smiles:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                ik = Chem.MolToInchiKey(mol)
                if ik:
                    ik_full[ik] = cid
                    ik_b1[ik.split("-")[0]].append(cid)
    cw_db, cw_pc = {}, {}
    for cid, db_id, pc_cid in conn_local.execute("SELECT chembl_id, drugbank_id, pubchem_cid FROM crosswalk").fetchall():
        if db_id and db_id.strip(): cw_db[db_id.strip()] = cid
        if pc_cid and str(pc_cid).strip(): cw_pc[str(pc_cid).strip()] = cid
    conn_local.close()

    # DrugCentral structs
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

    # Rejected sids
    rejected_sids = set()
    with open(DATA_DIR / "rejected_rule_b.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rejected_sids.add(int(r["dc_struct_id"]))

    def map_one(sid, ik):
        clean_ik = ik.strip() if ik else None
        b1 = clean_ik.split("-")[0] if clean_ik else None
        if clean_ik and clean_ik in ik_full:
            return ik_full[clean_ik], "a"
        elif clean_ik and b1 in ik_b1:
            return ik_b1[b1][0], "b"
        elif sid in dc_idents:
            ids = dc_idents[sid]
            for db_id in ids.get("DRUGBANK_ID", []):
                if db_id in cw_db: return cw_db[db_id], "c"
            for pc_id in ids.get("PUBCHEM_CID", []):
                if pc_id in cw_pc: return cw_pc[pc_id], "c"
            for ch_id in ids.get("ChEMBL_ID", []):
                if ch_id in drugs: return ch_id, "c"
        return None, None

    struct_to_cid = {}
    for sid, name, ik in dc_ind_structs:
        cid, rule = map_one(sid, ik)
        if cid:
            struct_to_cid[sid] = (cid, rule)

    excluded_ind_ids = {r["indication_id"] for r in csv.DictReader(open(DATA_DIR / "excluded_indications.csv", encoding="utf-8"))}

    raw_indications = []
    with open(DATA_DIR / "drugcentral_raw_approved_indications.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sid = int(r["struct_id"])
            if sid in struct_to_cid:
                cid, rule = struct_to_cid[sid]
                raw_indications.append((sid, cid, rule, r["concept_id"], r["concept_name"]))

    pool_primary = {cid for cid, d in drugs.items() if not d["triv"]} - INORGANIC_EXCLUDE
    pool_with_inorganics = {cid for cid, d in drugs.items() if not d["triv"]}

    # C3: Rule (a) only, primary pool
    pairs_c3 = set()
    for sid, cid, rule, ind_id, ind_name in raw_indications:
        if rule == "a" and cid in pool_primary and ind_id not in excluded_ind_ids:
            pairs_c3.add((cid, ind_id, ind_name))

    pairs_c3_sorted = sorted(list(pairs_c3))
    out_c3 = DATA_DIR / "gold_pairs_rule_a_only.csv"
    with open(out_c3, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["chembl_id", "indication_id", "indication_name"])
        w.writerows(pairs_c3_sorted)
    print(f"C3 Gold (rule a only): {len(pairs_c3_sorted)} pairs across {len({p[0] for p in pairs_c3_sorted})} drugs -> {out_c3.name}")

    # C4: DB-A with 10 inorganics, all accepted rules
    pairs_c4 = set()
    for sid, cid, rule, ind_id, ind_name in raw_indications:
        if sid not in rejected_sids and cid in pool_with_inorganics and ind_id not in excluded_ind_ids:
            pairs_c4.add((cid, ind_id, ind_name))

    pairs_c4_sorted = sorted(list(pairs_c4))
    out_c4 = DATA_DIR / "gold_pairs_with_inorganics.csv"
    with open(out_c4, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["chembl_id", "indication_id", "indication_name"])
        w.writerows(pairs_c4_sorted)
    print(f"C4 Gold (with inorganics): {len(pairs_c4_sorted)} pairs across {len({p[0] for p in pairs_c4_sorted})} drugs -> {out_c4.name}")

if __name__ == "__main__":
    main()

