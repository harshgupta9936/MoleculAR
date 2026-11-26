import sqlite3

def main():
    print("Updating indications in drugrepo.db...")
    repo_conn = sqlite3.connect('drugrepo.db')
    repo_c = repo_conn.cursor()
    
    chembl_db_path = os.getenv("CHEMBL_DB_PATH", os.path.join(os.path.dirname(__file__), "..", "chembl_37", "chembl_37_sqlite", "chembl_37.db"))
    chembl_conn = sqlite3.connect(chembl_db_path)
    chembl_c = chembl_conn.cursor()
    
    # Get all drugs that have 'Other' or empty indication
    repo_c.execute("SELECT chembl_id, generic_base_name FROM drugs WHERE indication = 'Other' OR indication IS NULL OR indication = ''")
    drugs_to_update = repo_c.fetchall()
    
    updates = []
    found = 0
    
    for cid, name in drugs_to_update:
        if str(cid).startswith('CHEMBL'):
            chembl_c.execute('''
                SELECT i.mesh_heading 
                FROM drug_indication i 
                JOIN molecule_dictionary m ON i.molregno = m.molregno 
                WHERE m.chembl_id = ?
            ''', (cid,))
            rows = chembl_c.fetchall()
            if rows:
                inds = set([r[0].lower() for r in rows])
                ind_str = ", ".join(sorted(inds))
                updates.append((ind_str, cid))
                found += 1
                
    print(f"Found real indications for {found} out of {len(drugs_to_update)} drugs from ChEMBL.")
    
    repo_c.executemany("UPDATE drugs SET indication = ? WHERE chembl_id = ?", updates)
    repo_conn.commit()
    
    # Check how many are still 'Other'
    repo_c.execute("SELECT COUNT(*) FROM drugs WHERE indication = 'Other'")
    print(f"Drugs still labeled 'Other': {repo_c.fetchone()[0]}")
    
    repo_conn.close()
    chembl_conn.close()
    print("Database updated.")

if __name__ == "__main__":
    main()
