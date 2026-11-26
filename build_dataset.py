import os
from chembl_webresource_client.new_client import new_client
import pandas as pd
import json
import time
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

salt_remover = SaltRemover()

def validate_and_canonicalize(smiles):
    """
    Validate a SMILES string, strip salts, and return the canonical form.
    Returns None if the SMILES is invalid or empty after processing.
    """
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Strip salt fragments (keep largest organic fragment)
    try:
        mol = salt_remover.StripMol(mol)
    except Exception:
        pass
    # Re-sanitize after stripping
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    canonical = Chem.MolToSmiles(mol)
    if not canonical:
        return None, None
    inchikey = Chem.MolToInchiKey(mol)
    return canonical, inchikey


# ATC code -> therapeutic class name mapping
ATC_TO_INDICATION = {
    "C09": "Cardiovascular / Hypertension",
    "A10": "Diabetes / Metabolic",
    "N06A": "Depression / Antidepressant",
    "M01A": "Pain / Anti-inflammatory",
    "C10": "Cardiovascular / Lipid-lowering",
    "N03A": "Epilepsy / Anticonvulsant",
    "G04C": "Urological / BPH",
    "N05A": "Psychosis / Antipsychotic",
}

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

    atc_codes = ["C09", "A10", "N06A", "M01A", "C10", "N03A", "G04C", "N05A"]
    records = []
    
    print("Fetching Phase 4 drugs and filtering by ATC codes... (This may take a few minutes)")
    approved_drugs = molecule.filter(max_phase=4).only(
        ["molecule_chembl_id", "pref_name", "molecule_structures", "atc_classifications"]
    )
    
    # Track how many we have found for each code
    found_counts = {code: 0 for code in atc_codes}
    smiles_dropped = 0
    
    for r in approved_drugs:
        if r.get("molecule_structures") and r.get("pref_name") and r.get("atc_classifications"):
            raw_smiles = r["molecule_structures"].get("canonical_smiles")
            res = validate_and_canonicalize(raw_smiles)
            if res == (None, None) or res is None:
                smiles_dropped += 1
                continue
            canonical, inchikey = res
            
            for code in atc_codes:
                if any(atc.startswith(code) for atc in r["atc_classifications"]):
                    drug_name = r["pref_name"].upper()
                    records.append({
                        "drug_name": drug_name,
                        "chembl_id": r["molecule_chembl_id"],
                        "smiles": canonical,
                        "inchikey": inchikey,
                        "atc_category": code,
                        "indication": ATC_TO_INDICATION.get(code, "Other"),
                    })
                    break # Drug assigned to this code

    print(f"SMILES validation: {smiles_dropped} entries dropped due to invalid/unparseable structures")

    # Must-have drugs for validation cases
    validation_drugs = ["SILDENAFIL", "MINOXIDIL", "METFORMIN", "THALIDOMIDE", 
                        "DULOXETINE", "GABAPENTIN", "BUPROPION", "RALOXIFENE"]

    existing_names = {r["drug_name"] for r in records}
    missing_drugs = [d for d in validation_drugs if d not in existing_names]

    if missing_drugs:
        print(f"\nFetching explicitly requested validation drugs: {missing_drugs}")
        for drug in missing_drugs:
            res = molecule.filter(pref_name__iexact=drug).only(
                ["molecule_chembl_id", "pref_name", "molecule_structures"]
            )
            for r in res:
                if r.get("molecule_structures") and r.get("pref_name"):
                    name = r["pref_name"].upper()
                    if name not in existing_names:
                        res = validate_and_canonicalize(raw_smiles)
                        if res == (None, None) or res is None:
                            continue
                        canonical, inchikey = res
                        records.append({
                            "drug_name": name,
                            "chembl_id": r["molecule_chembl_id"],
                            "smiles": canonical,
                            "inchikey": inchikey,
                            "atc_category": "VALIDATION",
                            "indication": VALIDATION_INDICATIONS.get(name, "Unknown"),
                        })
                        existing_names.add(name)
                        break

    # For drugs already in the dataset that are also validation drugs, patch their indication
    for rec in records:
        if rec["drug_name"] in VALIDATION_INDICATIONS:
            rec["indication"] = VALIDATION_INDICATIONS[rec["drug_name"]]

    df_new = pd.DataFrame(records)
    
    # Track raw count before deduplication
    raw_count = len(df_new)
    print(f"\nRaw pull count (Phase 4, matching ATC codes): {raw_count}")
    
    # Name deduplication just in case
    df_new = df_new.drop_duplicates(subset="drug_name")

    # Final deduplication by InChIKey (catches same molecule under different names/salts)
    before_dedup = len(df_new)
    df_new = df_new.drop_duplicates(subset="inchikey", keep="first")
    after_dedup = len(df_new)
    if raw_count != after_dedup:
        print(f"Deduplication: removed {raw_count - after_dedup} duplicate structures/names")

    print(f"\nTotal unique drugs in expanded dataset (post-dedup): {len(df_new)}")
    df_new.to_csv("full_dataset.csv", index=False)
    print("Saved to full_dataset.csv")

    # --- Fetch targets with human-readable names ---
    print("\nFetching target mechanisms for drugs...")
    drug_targets = {}
    
    # Collect all unique target IDs first for batch resolution
    all_target_ids = set()

    for i, row in df_new.iterrows():
        chembl_id = row["chembl_id"]
        time.sleep(0.1)
        try:
            mech_res = mechanism.filter(molecule_chembl_id=chembl_id)
            targets = [m.get("target_chembl_id") for m in mech_res if m.get("target_chembl_id")]
            unique_targets = list(set(targets))
            drug_targets[chembl_id] = unique_targets
            all_target_ids.update(unique_targets)
        except Exception as e:
            print(f"Error fetching targets for {chembl_id}: {e}")
            drug_targets[chembl_id] = []

    # Resolve target IDs to human-readable names
    print(f"\nResolving {len(all_target_ids)} unique target IDs to names...")
    target_name_map = {}
    for tid in all_target_ids:
        time.sleep(0.05)
        try:
            t = target_api.get(tid)
            if t and t.get("pref_name"):
                target_name_map[tid] = t["pref_name"]
            else:
                target_name_map[tid] = tid  # Fallback to ID
        except Exception:
            target_name_map[tid] = tid

    # Build enriched targets: {"CHEMBL_ID": [{"target_id": "...", "target_name": "..."}, ...]}
    drug_targets_enriched = {}
    for chembl_id, tids in drug_targets.items():
        drug_targets_enriched[chembl_id] = [
            {"target_id": tid, "target_name": target_name_map.get(tid, tid)}
            for tid in tids
        ]

    with open("drug_targets.json", "w") as f:
        json.dump(drug_targets_enriched, f, indent=2)

    # Also save the flat target name map for quick lookups
    with open("target_name_map.json", "w") as f:
        json.dump(target_name_map, f, indent=2)

    drugs_with_targets = sum(1 for v in drug_targets.values() if v)
    target_coverage = drugs_with_targets / len(df_new)
    print(f"\nTarget coverage: {target_coverage*100:.1f}% ({drugs_with_targets}/{len(df_new)} drugs have known targets)")
    print("Saved drug_targets.json and target_name_map.json")

if __name__ == "__main__":
    main()
