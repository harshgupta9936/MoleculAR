import sqlite3
import pandas as pd
import json
import time

def build_database():
    print("Initializing SQLite database...")
    conn = sqlite3.connect("drugrepo.db")
    cursor = conn.cursor()
    
    # Create tables
    cursor.execute("DROP TABLE IF EXISTS drugs")
    cursor.execute("DROP TABLE IF EXISTS targets")
    cursor.execute("DROP TABLE IF EXISTS side_effects")
    
    cursor.execute("""
        CREATE TABLE drugs (
            chembl_id TEXT PRIMARY KEY,
            drug_name TEXT,
            smiles TEXT,
            atc_category TEXT,
            indication TEXT,
            fingerprint_blob TEXT,
            hubness_score REAL,
            hubness_percentile REAL
        )
    """)
    
    cursor.execute("""
        CREATE TABLE targets (
            chembl_id TEXT,
            uniprot_id TEXT,
            target_name TEXT,
            FOREIGN KEY(chembl_id) REFERENCES drugs(chembl_id)
        )
    """)
    
    cursor.execute("""
        CREATE TABLE side_effects (
            chembl_id TEXT,
            side_effect_name TEXT,
            FOREIGN KEY(chembl_id) REFERENCES drugs(chembl_id)
        )
    """)
    
    # Indexes
    cursor.execute("CREATE INDEX idx_drug_name ON drugs(drug_name)")
    
    print("Loading data from colab_output_dataset.csv...")
    try:
        df = pd.read_csv("colab_output_dataset.csv", dtype={'morgan_fp': str})
    except FileNotFoundError:
        print("colab_output_dataset.csv not found! Run the colab_rebuild.py in Colab first.")
        return
        
    print(f"Loaded {len(df)} deduplicated, fingerprinted drugs.")
    
    # Insert drugs
    drugs_to_insert = []
    for _, row in df.iterrows():
        drugs_to_insert.append((
            row["chembl_id"],
            row["drug_name"],
            row["desalted_smiles"],
            row["atc_category"],
            row["indication"],
            row["morgan_fp"],
            row["hubness_score"],
            row["hubness_percentile"]
        ))
        
    cursor.executemany(
        "INSERT INTO drugs (chembl_id, drug_name, smiles, atc_category, indication, fingerprint_blob, hubness_score, hubness_percentile) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        drugs_to_insert
    )
    
    print("Loading targets from drug_targets_robust.json...")
    try:
        with open("drug_targets_robust.json", "r") as f:
            targets_data = json.load(f)
    except FileNotFoundError:
        print("drug_targets_robust.json not found!")
        targets_data = {}

    BLACKLISTED_TARGETS = [
        "LMNA", "Prelamin A/C", # Artifact targets
        "Titin", "TTN",
        "Zinc finger",
        "MUC", "Mucin"
    ]
        
    targets_inserted = 0
    targets_dropped = 0
    for chembl_id, t_info in targets_data.items():
        if chembl_id not in df['chembl_id'].values:
            continue
            
        targets = t_info.get("targets", [])
        for t in targets:
            name = t.get("target_name", "")
            if any(b.lower() in name.lower() for b in BLACKLISTED_TARGETS):
                targets_dropped += 1
                continue
                
            cursor.execute(
                "INSERT INTO targets (chembl_id, uniprot_id, target_name) VALUES (?, ?, ?)",
                (chembl_id, t.get("uniprot_id", ""), name)
            )
            targets_inserted += 1

    print(f"Inserted {targets_inserted} targets. Dropped {targets_dropped} artifact targets.")
    
    print("Loading side effects from drug_side_effects.json...")
    try:
        with open("drug_side_effects.json", "r") as f:
            se_data = json.load(f)
    except FileNotFoundError:
        print("drug_side_effects.json not found!")
        se_data = {}
        
    SUFFIXES = [
        " HYDROCHLORIDE", " SODIUM", " POTASSIUM", " CALCIUM", " MAGNESIUM", " ZINC",
        " SULFATE", " MESYLATE", " BESYLATE", " TOSYLATE", " MALEATE", " FUMARATE",
        " TARTRATE", " CITRATE", " ACETATE", " PHOSPHATE", " DIPHOSPHATE", " NITRATE",
        " BROMIDE", " HYDROBROMIDE", " CHLORIDE", " IODIDE", " FLUORIDE", " LACTATE", " SUCCINATE",
        " PAMOATE", " VALERATE", " STEARATE", " PALMITATE", " OLEATE", " SALICYLATE",
        " BENZOATE", " NAPSYLATE", " GLUCONATE", " GLUCEPTATE", " LACTOBIONATE",
        " ANHYDROUS", " MONOHYDRATE", " DIHYDRATE", " TRIHYDRATE", " TETRAHYDRATE",
        " PENTA", " HEXA", " HEPTA", " OCTA", " NAI", " FREE ACID", " FREE BASE",
        " CATION", " ANION", ", (E)-", ", (Z)-"
    ]

    def clean_drug_name(name):
        cleaned = str(name).upper()
        changed = True
        while changed:
            changed = False
            for suffix in SUFFIXES:
                if cleaned.endswith(suffix):
                    cleaned = cleaned[:-len(suffix)].strip()
                    changed = True
        return cleaned

    name_to_cid = {clean_drug_name(row["drug_name"]): row["chembl_id"] for _, row in df.iterrows()}
    se_inserted = 0
    for d_name, se_info in se_data.items():
        clean_name = str(d_name).upper()
        if clean_name not in name_to_cid:
            continue
            
        cid = name_to_cid[clean_name]
        ses = se_info
        for se in ses:
            cursor.execute(
                "INSERT INTO side_effects (chembl_id, side_effect_name) VALUES (?, ?)",
                (cid, se)
            )
            se_inserted += 1
            
    print(f"Inserted {se_inserted} side effects.")
    
    conn.commit()
    conn.close()
    print("Successfully built drugrepo.db!")

if __name__ == "__main__":
    build_database()
