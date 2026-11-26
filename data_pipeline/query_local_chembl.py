import sqlite3
import pandas as pd
import sys

ATC_TO_INDICATION = {
    "C09": "Cardiovascular / Hypertension",
    "A10": "Diabetes / Metabolic",
    "N06A": "Depression / Antidepressant",
    "M01A": "Pain / Anti-inflammatory",
    "C10": "Cardiovascular / Lipid-lowering",
    "N03A": "Epilepsy / Anticonvulsant",
    "G04C": "Urological / BPH",
    "N05A": "Psychosis / Antipsychotic",
}

VALIDATION_INDICATIONS = {
    "SILDENAFIL": "Erectile dysfunction (originally Angina)",
    "MINOXIDIL": "Hair loss (originally Hypertension)",
    "METFORMIN": "Type 2 Diabetes",
    "THALIDOMIDE": "Multiple myeloma / Leprosy (originally Sedative)",
    "DULOXETINE": "Depression / Diabetic neuropathy",
    "GABAPENTIN": "Epilepsy / Neuropathic pain",
    "BUPROPION": "Depression / Smoking cessation",
    "RALOXIFENE": "Osteoporosis / Breast cancer prevention",
}

def build_dataset_from_sqlite(db_path):
    print(f"Connecting to {db_path}...")
    conn = sqlite3.connect(db_path)
    
    query = """
    SELECT md.chembl_id, md.pref_name, cs.canonical_smiles, cs.standard_inchi_key, 
           GROUP_CONCAT(mac.level5) as atc_codes
    FROM molecule_dictionary md
    JOIN compound_structures cs ON md.molregno = cs.molregno
    LEFT JOIN molecule_atc_classification mac ON md.molregno = mac.molregno
    WHERE md.max_phase = 4
    GROUP BY md.chembl_id, md.pref_name, cs.canonical_smiles, cs.standard_inchi_key
    """
    
    print("Executing query for all max_phase=4 drugs...")
    df_raw = pd.read_sql_query(query, conn)
    print(f"Raw rows from DB: {len(df_raw)}")
    
    records = []
    atc_codes = ["C09", "A10", "N06A", "M01A", "C10", "N03A", "G04C", "N05A"]
    
    seen_chembl = set()
    for _, row in df_raw.iterrows():
        cid = row['chembl_id']
        if cid in seen_chembl: continue
        
        raw_smiles = row['canonical_smiles']
        pref_name = row['pref_name']
        atc_str = row['atc_codes']
        inchikey = row['standard_inchi_key']
        
        if not pref_name or not raw_smiles or not inchikey:
            continue
            
        atc_list = atc_str.split(',') if pd.notnull(atc_str) else []
        assigned_code = "Other"
        for code in atc_codes:
            if any(atc.startswith(code) for atc in atc_list):
                assigned_code = code
                break
                
        records.append({
            "drug_name": pref_name.upper(),
            "chembl_id": cid,
            "smiles": raw_smiles,
            "inchikey": inchikey,
            "atc_category": assigned_code,
            "indication": ATC_TO_INDICATION.get(assigned_code, "Other"),
        })
        seen_chembl.add(cid)
        
    validation_drugs = list(VALIDATION_INDICATIONS.keys())
    existing_names = {r["drug_name"] for r in records}
    missing_drugs = [d for d in validation_drugs if d not in existing_names]
    
    if missing_drugs:
        print(f"\nFetching explicitly requested validation drugs: {missing_drugs}")
        placeholders = ','.join(['?']*len(missing_drugs))
        v_query = f"""
        SELECT md.chembl_id, md.pref_name, cs.canonical_smiles, cs.standard_inchi_key
        FROM molecule_dictionary md
        JOIN compound_structures cs ON md.molregno = cs.molregno
        WHERE UPPER(md.pref_name) IN ({placeholders})
        """
        df_val = pd.read_sql_query(v_query, conn, params=missing_drugs)
        for _, row in df_val.iterrows():
            cid = row['chembl_id']
            pref_name = row['pref_name'].upper()
            raw_smiles = row['canonical_smiles']
            inchikey = row['standard_inchi_key']
            
            if not raw_smiles or not inchikey: continue
                
            records.append({
                "drug_name": pref_name,
                "chembl_id": cid,
                "smiles": raw_smiles,
                "inchikey": inchikey,
                "atc_category": "VALIDATION",
                "indication": VALIDATION_INDICATIONS.get(pref_name, "Unknown"),
            })
            existing_names.add(pref_name)

    conn.close()

    for rec in records:
        if rec["drug_name"] in VALIDATION_INDICATIONS:
            rec["indication"] = VALIDATION_INDICATIONS[rec["drug_name"]]
            
    df_new = pd.DataFrame(records)
    raw_count = len(df_new)
    print(f"\nRaw pull count (Phase 4, matching ATC codes + validation): {raw_count}")
    
    df_new = df_new.drop_duplicates(subset="inchikey", keep="first")
    after_dedup = len(df_new)
    if raw_count != after_dedup:
        print(f"Deduplication: removed {raw_count - after_dedup} duplicate structures")

    print(f"\nTotal unique drugs in expanded dataset (post-dedup): {after_dedup}")
    df_new.to_csv("full_dataset.csv", index=False)
    print("Saved to full_dataset.csv")

    benchmarks = ['SILDENAFIL', 'MINOXIDIL', 'METFORMIN', 'THALIDOMIDE', 'DULOXETINE', 'GABAPENTIN', 'BUPROPION', 'RALOXIFENE']
    missing = [d for d in benchmarks if not df_new['drug_name'].str.contains(d).any()]
    if not missing:
        print("\nSUCCESS: All 8 locked benchmark drugs are present in the rebuilt dataset.")
    else:
        print(f"\nWARNING: Missing benchmark drugs: {missing}")

if __name__ == "__main__":
    build_dataset_from_sqlite("chembl_37/chembl_37_sqlite/chembl_37.db")
