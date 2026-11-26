import os
import sys
from chembl_webresource_client.new_client import new_client
import pandas as pd
import json
import time
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

salt_remover = SaltRemover()

from rdkit.Chem.MolStandardize import rdMolStandardize

largest_fragment = rdMolStandardize.LargestFragmentChooser()

def validate_and_canonicalize(smiles):
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        # Extract the largest covalently bonded fragment (drops salts/solvates)
        mol = largest_fragment.choose(mol)
    except Exception:
        pass
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    canonical = Chem.MolToSmiles(mol)
    if not canonical:
        return None
    return canonical

# Hardcoded indications for the 8 validation drugs (known ground truth)
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

def main():
    molecule = new_client.molecule
    mechanism = new_client.mechanism
    target_api = new_client.target

    records = []
    
    print("Fetching ALL approved drugs (max_phase=4)... (This will take a moment)")
    approved_drugs = molecule.filter(max_phase=4).only(
        ["molecule_chembl_id", "pref_name", "molecule_structures"]
    )
    
    smiles_dropped = 0
    total_pulled = 0
    
    for r in approved_drugs:
        total_pulled += 1
        if total_pulled % 500 == 0:
            print(f"Processed {total_pulled} drug records...")
            
        if r.get("molecule_structures") and r.get("pref_name"):
            raw_smiles = r["molecule_structures"].get("canonical_smiles")
            canonical = validate_and_canonicalize(raw_smiles)
            if canonical is None:
                smiles_dropped += 1
                continue
            
            drug_name = r["pref_name"].upper()
            records.append({
                "drug_name": drug_name,
                "chembl_id": r["molecule_chembl_id"],
                "smiles": canonical,
                "indication": VALIDATION_INDICATIONS.get(drug_name, "Approved Drug"),
            })

    print(f"Total raw molecules fetched: {total_pulled}")
    print(f"SMILES validation: {smiles_dropped} entries dropped due to invalid/unparseable structures")

    df_new = pd.DataFrame(records).drop_duplicates(subset="drug_name")
    
    # Final deduplication by structural InChIKey (first 25 chars for skeleton + stereochem)
    before_dedup = len(df_new)
    
    def get_inchi_key(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        inchi = Chem.MolToInchiKey(mol)
        if not inchi: return None
        return inchi[:25]  # Drop only the protonation state character (last char) and version
        
    df_new["inchikey_25"] = df_new["smiles"].apply(get_inchi_key)
    df_new = df_new.drop_duplicates(subset="inchikey_25", keep="first")
    df_new = df_new.drop(columns=["inchikey_25"])
    
    after_dedup = len(df_new)
    if before_dedup != after_dedup:
        print(f"Deduplication: removed {before_dedup - after_dedup} structural duplicates (salts/isomers/tautomers)")

    print(f"\nTotal unique drugs in expanded dataset: {len(df_new)}")
    df_new.to_csv("full_dataset.csv", index=False)
    print("Saved to full_dataset.csv")

    # --- Fetch targets ---
    print("\nFetching target mechanisms for drugs...")
    drug_targets = {}
    all_target_ids = set()

    total_drugs = len(df_new)
    for i, row in df_new.iterrows():
        chembl_id = row["chembl_id"]
        if i % 100 == 0:
            print(f"Fetching targets: {i}/{total_drugs}")
            
        # Optional: no sleep needed usually for chembl_webresource_client caching, but adding tiny delay just in case
        time.sleep(0.02)
        try:
            # We filter for only the drug's mechanisms
            mech_res = mechanism.filter(molecule_chembl_id=chembl_id)
            targets = []
            for m in mech_res:
                # Keep only high confidence (if available) or direct targets
                # 'mechanism_of_action' endpoint provides direct targets.
                if m.get("target_chembl_id"):
                    targets.append(m.get("target_chembl_id"))
            
            unique_targets = list(set(targets))
            drug_targets[chembl_id] = unique_targets
            all_target_ids.update(unique_targets)
        except Exception as e:
            # Silent continue for network blips
            drug_targets[chembl_id] = []

    # Resolve target IDs to human-readable names
    print(f"\nResolving {len(all_target_ids)} unique target IDs to names...")
    target_name_map = {}
    for idx, tid in enumerate(all_target_ids):
        if idx % 50 == 0:
            print(f"Resolving targets: {idx}/{len(all_target_ids)}")
        time.sleep(0.02)
        try:
            t = target_api.get(tid)
            if t and t.get("pref_name"):
                target_name_map[tid] = t["pref_name"]
            else:
                target_name_map[tid] = tid
        except Exception:
            target_name_map[tid] = tid

    drug_targets_enriched = {}
    for chembl_id, tids in drug_targets.items():
        drug_targets_enriched[chembl_id] = [
            {"target_id": tid, "target_name": target_name_map.get(tid, tid)}
            for tid in tids
        ]

    with open("drug_targets.json", "w") as f:
        json.dump(drug_targets_enriched, f, indent=2)

    with open("target_name_map.json", "w") as f:
        json.dump(target_name_map, f, indent=2)

    drugs_with_targets = sum(1 for v in drug_targets.values() if v)
    target_coverage = drugs_with_targets / len(df_new)
    print(f"\nTarget coverage: {target_coverage*100:.1f}% ({drugs_with_targets}/{len(df_new)} drugs have known targets)")
    print("Saved drug_targets.json and target_name_map.json")

if __name__ == "__main__":
    main()
