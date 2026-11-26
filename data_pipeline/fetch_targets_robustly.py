import json
import pandas as pd
import time
from chembl_webresource_client.new_client import new_client
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

def fetch_drug(row):
    cid = row["chembl_id"]
    drug_name = row["drug_name"]
    
    success = False
    mechs = []
    
    # 1. Retry mechanism first
    for attempt in range(3):
        try:
            url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id={cid}"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                mechs = data.get("mechanisms", [])
                success = True
                break
            else:
                time.sleep(2 ** attempt)
        except Exception:
            time.sleep(2 ** attempt)
            
    if success and len(mechs) > 0:
        targets = list(set([m.get("target_chembl_id") for m in mechs if m.get("target_chembl_id")]))
        return cid, {"source": "mechanism", "target_ids": targets}
    else:
        # 2. Fallback to activity
        try:
            url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={cid}&pchembl_value__gte=6.0"
            resp = requests.get(url, timeout=15)
            if resp.status_code == 200:
                data = resp.json()
                activities = data.get("activities", [])
                # Explicitly filter for Binding (B) or Functional (F) assays to exclude ADMET (A) and Toxicity (T) artifacts
                targets = list(set([a.get("target_chembl_id") for a in activities if a.get("target_chembl_id") and a.get("assay_type") in ["B", "F"]]))
                return cid, {"source": "activity_fallback", "target_ids": targets}
            else:
                return cid, {"source": "activity_fallback_failed", "target_ids": []}
        except Exception:
            return cid, {"source": "activity_fallback_failed", "target_ids": []}

def fetch_targets_robustly():
    df = pd.read_csv("full_dataset.csv")
    total = len(df)
    print(f"Fetching robust targets for {total} drugs using ThreadPool...")
    
    drug_targets = {}
    all_target_ids = set()
    
    with ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(fetch_drug, row): i for i, row in df.iterrows()}
        completed = 0
        for future in as_completed(futures):
            completed += 1
            if completed % 100 == 0:
                print(f"Processed {completed}/{total} drugs...")
            cid, data = future.result()
            drug_targets[cid] = data
            all_target_ids.update(data["target_ids"])
            
    # 3. Resolve UniProt
    print(f"\nResolving {len(all_target_ids)} unique target IDs to UniProt accessions...")
    target_uniprot_map = {}
    target_api = new_client.target
    
    for idx, tid in enumerate(list(all_target_ids)):
        if idx % 50 == 0:
            print(f"Resolving targets: {idx}/{len(all_target_ids)}")
        time.sleep(0.05)
        try:
            t = target_api.get(tid)
            components = t.get("target_components", [])
            accessions = []
            component_names = []
            if components:
                for comp in components:
                    acc = comp.get("accession")
                    if acc: accessions.append(acc)
                    c_desc = comp.get("component_description")
                    if c_desc: component_names.append(c_desc)
            
            pref = t.get("pref_name", tid)
            if not component_names:
                component_names = [pref]
                
            target_uniprot_map[tid] = {
                "uniprot_accessions": accessions,
                "name": " / ".join(component_names) if component_names else pref
            }
        except Exception:
            target_uniprot_map[tid] = {"uniprot_accessions": [], "name": tid}

    # Build final
    final_output = {}
    for cid, data in drug_targets.items():
        enriched_targets = []
        for tid in data["target_ids"]:
            mapped = target_uniprot_map.get(tid, {})
            if mapped.get("uniprot_accessions"):
                for acc in mapped["uniprot_accessions"]:
                    enriched_targets.append({
                        "uniprot_id": acc,
                        "target_name": mapped["name"]
                    })
            else:
                enriched_targets.append({
                    "target_chembl_id": tid,
                    "target_name": mapped.get("name", tid)
                })
        
        # Deduplicate
        seen = set()
        unique_enriched = []
        for et in enriched_targets:
            key = et.get("uniprot_id") or et.get("target_chembl_id")
            if key not in seen:
                seen.add(key)
                unique_enriched.append(et)
                
        final_output[cid] = {
            "source": data["source"],
            "targets": unique_enriched
        }

    with open("drug_targets_robust.json", "w") as f:
        json.dump(final_output, f, indent=2)

    mech_count = sum(1 for v in final_output.values() if v["source"] == "mechanism")
    act_count = sum(1 for v in final_output.values() if v["source"] == "activity_fallback")
    fail_count = sum(1 for v in final_output.values() if "failed" in v["source"])
    
    print("\n--- FINAL STATS ---")
    print(f"Total Drugs: {total}")
    print(f"Mechanism Source: {mech_count}")
    print(f"Activity Fallback Source: {act_count}")
    print(f"Failed completely: {fail_count}")
    print("Done!")

if __name__ == "__main__":
    fetch_targets_robustly()
