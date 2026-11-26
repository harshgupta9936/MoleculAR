import time
import requests
import subprocess
import sys
import pandas as pd

def check_api():
    try:
        resp = requests.get('https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4&limit=1', timeout=15)
        if resp.status_code == 200:
            # ensure it's not a 500 html page
            if 'html' not in resp.text[:100].lower():
                return True
    except Exception as e:
        print(f"API Check Exception: {e}")
    return False

def wait_for_api():
    backoffs = [300, 600, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900, 900] # ~3 hours
    for wait_time in backoffs:
        print(f"Checking ChEMBL API status...")
        if check_api():
            print("API is UP! Proceeding to rebuild.")
            return True
        print(f"API is down. Waiting {wait_time // 60} minutes before retrying...")
        time.sleep(wait_time)
        
    print("API is still down after maximum retries. Exiting.")
    return False

def check_benchmark():
    try:
        df = pd.read_csv('full_dataset.csv')
        benchmarks = ['SILDENAFIL', 'MINOXIDIL', 'METFORMIN', 'THALIDOMIDE', 'DULOXETINE', 'GABAPENTIN', 'BUPROPION', 'RALOXIFENE']
        missing = [d for d in benchmarks if not df['drug_name'].str.contains(d).any()]
        if not missing:
            print("\nSUCCESS: All 8 locked benchmark drugs are present in the rebuilt dataset.")
        else:
            print(f"\nWARNING: Missing benchmark drugs: {missing}")
    except Exception as e:
        print(f"Error checking benchmarks: {e}")

def main():
    if not wait_for_api():
        sys.exit(1)
        
    print("\nRunning build_dataset.py (InChIKey dedup, strict ATC 8 filtering)...")
    result = subprocess.run([sys.executable, "-u", "build_dataset.py"], capture_output=True, text=True)
    
    print("\n--- BUILD_DATASET OUTPUT ---")
    print(result.stdout)
    
    if result.stderr:
        print("--- ERRORS ---")
        print(result.stderr)
        
    if result.returncode == 0:
        check_benchmark()
        print("\nRebuild complete. Check the output above for raw vs post-dedup counts.")
    else:
        print("\nRebuild failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
