import json
import pandas as pd
import requests

# Whitelist of essential drugs for Phase 1 Validation
whitelist = {
    "SILDENAFIL", "VARDENAFIL", "TADALAFIL", "AVANAFIL",
    "MINOXIDIL", "FINASTERIDE", "DUTASTERIDE",
    "METFORMIN", "PHENFORMIN",
    "THALIDOMIDE", "LENALIDOMIDE", "POMALIDOMIDE",
    "DULOXETINE", "VENLAFAXINE", "MILNACIPRAN", "DESVENLAFAXINE",
    "GABAPENTIN", "PREGABALIN",
    "BUPROPION", "VARENICLINE", "CYTISINE",
    "RALOXIFENE", "TAMOXIFEN", "TOREMIFENE", "BAZEDOXIFENE",
    "CLOZAPINE", "OLANZAPINE", "QUETIAPINE", "NEFAZODONE", "FENOFIBRATE", "DICLOFENAC",
    "OFLOXACIN", "LEVOFLOXACIN"
}

def fetch_drug(cid):
    mechs = []
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id={cid}"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            mechs = data.get("mechanisms", [])
    except Exception as e:
        print(f"Error fetching {cid}: {e}")
        
    if len(mechs) > 0:
        targets = list(set([m.get("target_chembl_id") for m in mechs if m.get("target_chembl_id")]))
        return {"source": "mechanism", "target_ids": targets}
    
    # Fallback to activity
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={cid}&pchembl_value__gte=6.0"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            activities = data.get("activities", [])
            targets = list(set([a.get("target_chembl_id") for a in activities if a.get("target_chembl_id") and a.get("assay_type") in ["B", "F"]]))
            return {"source": "activity_fallback", "target_ids": targets}
    except Exception:
        pass
        
    return {"source": "failed", "target_ids": []}

def main():
    df = pd.read_csv("full_dataset.csv")
    final_output = {}
    
    # Fetch targets for whitelisted drugs
    count = 0
    for i, row in df.iterrows():
        cid = row["chembl_id"]
        name = row["drug_name"]
        
        # We must check if any whitelist string is IN the drug name (e.g. VARDENAFIL HYDROCHLORIDE)
        is_whitelisted = any(w in name for w in whitelist)
        
        if is_whitelisted:
            print(f"Fetching targets for {name} ({cid})...")
            data = fetch_drug(cid)
            enriched = [{"target_chembl_id": tid, "target_name": tid, "uniprot_id": tid} for tid in data["target_ids"]]
            final_output[cid] = {
                "source": data["source"],
                "targets": enriched
            }
            count += 1
        else:
            final_output[cid] = {
                "source": "skipped",
                "targets": []
            }
            
    print(f"Fetched targets for {count} whitelisted drugs. Writing to drug_targets_robust.json...")
    with open("drug_targets_robust.json", "w") as f:
        json.dump(final_output, f, indent=2)
        
    print("Done!")

if __name__ == "__main__":
    main()
