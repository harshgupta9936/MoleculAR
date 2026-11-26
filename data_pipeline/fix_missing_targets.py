import sqlite3

BLACKLISTED_TARGETS = [
    "LMNA", "Prelamin A/C", # Artifact targets
    "Titin", "TTN",
    "Zinc finger",
    "MUC", "Mucin"
]

def fix_missing_targets():
    repo_conn = sqlite3.connect("drugrepo.db")
    repo_cursor = repo_conn.cursor()
    
    chembl_db_path = os.getenv("CHEMBL_DB_PATH", os.path.join(os.path.dirname(__file__), "..", "chembl_37", "chembl_37_sqlite", "chembl_37.db"))
    chembl_conn = sqlite3.connect(chembl_db_path)
    chembl_cursor = chembl_conn.cursor()
    
    # 1. Find drugs with 0 targets
    repo_cursor.execute('''
    SELECT d.chembl_id FROM drugs d
    LEFT JOIN targets t ON d.chembl_id = t.chembl_id
    WHERE t.chembl_id IS NULL
    ''')
    zero_target_drugs = [row[0] for row in repo_cursor.fetchall()]
    print(f"Found {len(zero_target_drugs)} drugs with 0 targets.")
    
    targets_inserted = 0
    drugs_fixed = 0
    
    for cid in zero_target_drugs:
        # Get molregno
        chembl_cursor.execute("SELECT molregno FROM molecule_dictionary WHERE chembl_id=?", (cid,))
        row = chembl_cursor.fetchone()
        if not row:
            continue
        molregno = row[0]
        
        # Extract targets using hierarchy to catch salts/parents
        query = '''
        SELECT t.pref_name, c.accession
        FROM drug_mechanism dm
        JOIN target_dictionary t ON dm.tid = t.tid
        LEFT JOIN target_components tc ON t.tid = tc.tid
        LEFT JOIN component_sequences c ON tc.component_id = c.component_id
        WHERE dm.molregno IN (
            SELECT molregno FROM molecule_hierarchy 
            WHERE parent_molregno = (
                SELECT parent_molregno FROM molecule_hierarchy WHERE molregno = ?
            )
        )
        '''
        chembl_cursor.execute(query, (molregno,))
        targets = chembl_cursor.fetchall()
        
        if targets:
            drugs_fixed += 1
            # Insert unique targets
            unique_targets = {}
            for t_name, u_id in targets:
                if not t_name: continue
                # Skip blacklisted
                if any(b.lower() in t_name.lower() for b in BLACKLISTED_TARGETS):
                    continue
                unique_targets[t_name] = u_id
                
            for t_name, u_id in unique_targets.items():
                repo_cursor.execute("INSERT INTO targets (chembl_id, target_name, uniprot_id) VALUES (?, ?, ?)", (cid, t_name, u_id))
                targets_inserted += 1

    repo_conn.commit()
    repo_conn.close()
    chembl_conn.close()
    
    print(f"Fixed {drugs_fixed} drugs by inserting {targets_inserted} targets.")

if __name__ == "__main__":
    fix_missing_targets()
