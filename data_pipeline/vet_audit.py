import sqlite3
import requests
import time
import urllib.parse
import concurrent.futures
import threading

# Comprehensive list of salt and form suffixes to strip
SUFFIXES = [
    " HYDROCHLORIDE", " SODIUM", " POTASSIUM", " CALCIUM", " MAGNESIUM", " ZINC",
    " SULFATE", " MESYLATE", " BESYLATE", " TOSYLATE", " MALEATE", " FUMARATE",
    " TARTRATE", " CITRATE", " ACETATE", " PHOSPHATE", " DIPHOSPHATE", " NITRATE",
    " BROMIDE", " CHLORIDE", " IODIDE", " FLUORIDE", " LACTATE", " SUCCINATE",
    " PAMOATE", " VALERATE", " STEARATE", " PALMITATE", " OLEATE", " SALICYLATE",
    " BENZOATE", " NAPSYLATE", " GLUCONATE", " GLUCEPTATE", " LACTOBIONATE",
    " ANHYDROUS", " MONOHYDRATE", " DIHYDRATE", " TRIHYDRATE", " TETRAHYDRATE",
    " PENTA", " HEXA", " HEPTA", " OCTA", " NAI", " FREE ACID", " FREE BASE",
    " CATION", " ANION", ", (E)-", ", (Z)-"
]

def clean_drug_name(name):
    # Iteratively strip suffixes from the end until clean
    cleaned = name.upper()
    changed = True
    while changed:
        changed = False
        for suffix in SUFFIXES:
            if cleaned.endswith(suffix):
                cleaned = cleaned[:-len(suffix)].strip()
                changed = True
    return cleaned

def check_openfda(drug_name):
    base_name = clean_drug_name(drug_name)
    query = f'openfda.generic_name:"{base_name}"'
    url = f"https://api.fda.gov/drug/label.json?search={urllib.parse.quote(query)}&limit=1"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            return True
        elif r.status_code == 404:
            return False
    except requests.RequestException:
        pass
    return None

def check_dailymed(drug_name):
    base_name = clean_drug_name(drug_name)
    url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name={urllib.parse.quote(base_name)}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            if data.get("data") and len(data["data"]) > 0:
                return True
            else:
                return False
    except requests.RequestException:
        pass
    return None

def process_openfda_batch(chunk):
    results = []
    for cid, name in chunk:
        found = check_openfda(name)
        if found is False or found is None:
            results.append((cid, name))
        time.sleep(0.1) # Thread-local slight delay
    return results

def process_dailymed_batch(chunk):
    results = []
    for cid, name in chunk:
        found = check_dailymed(name)
        if not found:
            results.append((cid, name))
        time.sleep(0.1)
    return results

def chunk_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def main():
    conn = sqlite3.connect("drugrepo.db")
    cursor = conn.cursor()
    
    cursor.execute("SELECT chembl_id, drug_name FROM drugs")
    drugs = cursor.fetchall()
    
    print(f"Starting Veterinary/Human-use Audit for {len(drugs)} drugs (with concurrent threading & suffix stripping)...")
    
    openfda_misses = []
    
    # Step A: openFDA Triage (Multithreaded)
    print("Step A: Querying openFDA...")
    chunks = list(chunk_list(drugs, 100))
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(process_openfda_batch, chunk): chunk for chunk in chunks}
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            openfda_misses.extend(future.result())
            completed += len(futures[future])
            print(f"  Processed {completed}/{len(drugs)}...")
            
    print(f"openFDA found no human label for {len(openfda_misses)} drugs.")
    
    # Step B: DailyMed Verification (Multithreaded)
    print("Step B: Querying DailyMed for the zero-hits...")
    true_zero_hits = []
    
    dm_chunks = list(chunk_list(openfda_misses, 50))
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(process_dailymed_batch, chunk): chunk for chunk in dm_chunks}
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            true_zero_hits.extend(future.result())
            completed += len(futures[future])
            print(f"  Processed {completed}/{len(openfda_misses)}...")
        
    print(f"\nAudit complete. Found {len(true_zero_hits)} TRUE zero-hit compounds (missing in both APIs).")
    
    if true_zero_hits:
        print(f"{'ChEMBL ID':<15} | {'Drug Name'}")
        print("-" * 50)
        for cid, name in true_zero_hits:
            print(f"{cid:<15} | {name}")
    else:
        print("All drugs have human use labels confirmed.")
        
    print("\nPaused for review. The database has NOT been modified.")
    conn.close()

if __name__ == "__main__":
    main()
