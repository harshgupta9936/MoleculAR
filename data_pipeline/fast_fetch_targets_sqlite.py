import sqlite3
import pandas as pd
import json
import time

def fetch_targets_from_sqlite(db_path="chembl_37/chembl_37_sqlite/chembl_37.db"):
    print(f"Connecting to {db_path}...")
    conn = sqlite3.connect(db_path)
    
    # We will build a mapping of drug chembl_id -> list of targets
    # Query: for each drug in drug_mechanism, get the target info and uniprot accession
    
    query = """
    SELECT md.chembl_id as drug_chembl_id, 
           td.chembl_id as target_chembl_id, 
           td.pref_name as target_name, 
           cs.accession as uniprot_id
    FROM molecule_dictionary md
    JOIN drug_mechanism dm ON md.molregno = dm.molregno
    JOIN target_dictionary td ON dm.tid = td.tid
    LEFT JOIN target_components tc ON td.tid = tc.tid
    LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
    """
    print("Executing query for all drug mechanisms...")
    start = time.time()
    df_mech = pd.read_sql_query(query, conn)
    print(f"Query completed in {time.time() - start:.2f}s. Rows: {len(df_mech)}")
    conn.close()
    
    # Read the dataset we just built
    print("Reading full_dataset.csv...")
    try:
        df_dataset = pd.read_csv("full_dataset.csv")
    except FileNotFoundError:
        print("full_dataset.csv not found.")
        return
        
    final_output = {}
    
    # Pre-group mechanisms by drug_chembl_id for fast lookup
    print("Grouping mechanisms...")
    mech_groups = df_mech.groupby("drug_chembl_id")
    
    count = 0
    for cid in df_dataset["chembl_id"]:
        targets = []
        if cid in mech_groups.groups:
            group = mech_groups.get_group(cid)
            seen_keys = set()
            for _, row in group.iterrows():
                tid = row["target_chembl_id"]
                uid = row["uniprot_id"]
                tname = row["target_name"]
                
                # Key deduplication
                key = uid if pd.notnull(uid) else tid
                if key not in seen_keys:
                    seen_keys.add(key)
                    t_dict = {
                        "target_chembl_id": tid,
                        "target_name": tname
                    }
                    if pd.notnull(uid):
                        t_dict["uniprot_id"] = uid
                    targets.append(t_dict)
                    
        source = "mechanism" if len(targets) > 0 else "sqlite_fallback_failed"
        final_output[cid] = {
            "source": source,
            "targets": targets
        }
        if len(targets) > 0: count += 1

    print(f"Fetched targets for {count}/{len(df_dataset)} drugs in dataset. Writing to drug_targets_robust.json...")
    with open("drug_targets_robust.json", "w") as f:
        json.dump(final_output, f, indent=2)
        
    print("Done!")

if __name__ == "__main__":
    fetch_targets_from_sqlite()
