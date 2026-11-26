import sqlite3
import requests
import urllib.parse
import concurrent.futures
import time

# Suffixes to strip
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

def check_openfda(base_name):
    url = f'https://api.fda.gov/drug/label.json?search=openfda.generic_name:"{urllib.parse.quote(base_name)}"&limit=1'
    try:
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            return True
    except:
        pass
    return False

def check_dailymed(base_name):
    url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name={urllib.parse.quote(base_name)}"
    try:
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            data = r.json()
            if data.get("data") and len(data["data"]) > 0:
                return True
    except:
        pass
    return False

def check_synonyms(synonyms):
    cleaned_syns = list(set(clean_drug_name(s) for s in synonyms if s))
    for syn in cleaned_syns:
        if len(syn) < 3:
            continue
        if check_openfda(syn) or check_dailymed(syn):
            return syn
        time.sleep(0.05)
    return None

def process_drug(item):
    cid, name, veterinary, withdrawn, synonyms = item
    if cid == 'CHEMBL1316': # Carprofen manual override
        return (cid, name, 'Manual Review - Veterinary (Non-US)', None)
        
    rescued_name = check_synonyms(synonyms)
    if rescued_name:
        return (cid, name, 'Naming Mismatch', rescued_name)
    elif veterinary == 1:
        return (cid, name, 'Veterinary', None)
    elif withdrawn == 1:
        return (cid, name, 'Withdrawn', None)
    else:
        return (cid, name, 'International/Non-US', None)

def main():
    print("Adding market_status column to drugrepo.db...")
    repo_conn = sqlite3.connect('drugrepo.db')
    repo_cursor = repo_conn.cursor()
    
    try:
        repo_cursor.execute("ALTER TABLE drugs ADD COLUMN market_status VARCHAR(50) DEFAULT 'Active (US)'")
    except sqlite3.OperationalError:
        pass # Column already exists
    
    repo_conn.commit()
    
    # 1. Parse the 720 list
    md_path = os.getenv("VET_AUDIT_MD_PATH", os.path.join(os.path.dirname(__file__), "..", "archive_data", "vet_audit_flagged_list.md"))
    chembl_ids = []
    if os.path.exists(md_path):
        with open(md_path, 'r', encoding='utf-8') as f:
            in_table = False
            for line in f:
                if '|' in line and 'ChEMBL ID' in line:
                    in_table = True
                    continue
                if in_table and '|' in line and '---' not in line:
                    parts = line.split('|')
                    cid = parts[1].strip()
                    if cid.startswith('CHEMBL'):
                        chembl_ids.append(cid)

    # 2. Query ChEMBL for metadata and synonyms
    db_path = os.getenv("CHEMBL_DB_PATH", os.path.join(os.path.dirname(__file__), "..", "chembl_37", "chembl_37_sqlite", "chembl_37.db"))
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    tasks = []
    for cid in chembl_ids:
        cursor.execute("SELECT molregno, pref_name, veterinary, withdrawn_flag FROM molecule_dictionary WHERE chembl_id=?", (cid,))
        row = cursor.fetchone()
        if not row:
            continue
        molregno, pref_name, vet, withdrawn = row
        
        cursor.execute("SELECT synonyms FROM molecule_synonyms WHERE molregno=?", (molregno,))
        syns = [r[0] for r in cursor.fetchall()]
        if pref_name:
            syns.append(pref_name)
            
        tasks.append((cid, pref_name, vet, withdrawn, syns))
    
    conn.close()
    
    print(f"Categorizing {len(tasks)} drugs...")
    
    results = {'Naming Mismatch': [], 'Veterinary': [], 'Withdrawn': [], 'International/Non-US': [], 'Manual Review - Veterinary (Non-US)': []}
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        future_to_cid = {executor.submit(process_drug, item): item for item in tasks}
        completed = 0
        for future in concurrent.futures.as_completed(future_to_cid):
            res = future.result()
            cat = res[2]
            results[cat].append(res)
            completed += 1
            if completed % 50 == 0:
                print(f"  Processed {completed}/{len(tasks)}...")
                
    print("\nApplying updates to drugrepo.db...")
    for cat, items in results.items():
        print(f"{cat}: {len(items)}")
        for item in items:
            cid, name, category, rescued_name = item
            if category == 'Naming Mismatch':
                # Re-fold as Active (US) and correct generic_base_name
                repo_cursor.execute("UPDATE drugs SET market_status = 'Active (US)', generic_base_name = ? WHERE chembl_id = ?", (rescued_name, cid))
            else:
                repo_cursor.execute("UPDATE drugs SET market_status = ? WHERE chembl_id = ?", (category, cid))
                
    repo_conn.commit()
    repo_conn.close()
    print("Database migration complete.")

if __name__ == "__main__":
    main()
