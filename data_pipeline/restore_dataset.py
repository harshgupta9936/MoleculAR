import json
import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from chembl_webresource_client.new_client import new_client

molecule = new_client.molecule
largest_fragment = rdMolStandardize.LargestFragmentChooser()

print("Loading original 2,467 chembl_ids from drug_targets_robust.json...")
with open("drug_targets_robust.json", "r") as f:
    target_data = json.load(f)

cids = list(target_data.keys())
print(f"Found {len(cids)} drugs to process.")

def validate_and_canonicalize(smiles):
    if not smiles: return None
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    try: mol = largest_fragment.choose(mol)
    except: pass
    try: Chem.SanitizeMol(mol)
    except: return None
    return Chem.MolToSmiles(mol)

# Hardcoded indications for the 8 validation drugs
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

print("Fetching SMILES and names from ChEMBL...")
# We fetch in batches to avoid timeout
batch_size = 500
records = []
for i in range(0, len(cids), batch_size):
    batch_cids = cids[i:i+batch_size]
    print(f"  Fetching batch {i} to {i+len(batch_cids)}...")
    res = molecule.filter(molecule_chembl_id__in=batch_cids).only(["molecule_chembl_id", "pref_name", "molecule_structures"])
    for r in res:
        cid = r["molecule_chembl_id"]
        if r.get("molecule_structures") and r.get("pref_name"):
            raw_smiles = r["molecule_structures"].get("canonical_smiles")
            canonical = validate_and_canonicalize(raw_smiles)
            if not canonical: continue
            
            dname = r["pref_name"].upper()
            records.append({
                "drug_name": dname,
                "chembl_id": cid,
                "smiles": canonical,
                "indication": VALIDATION_INDICATIONS.get(dname, "Approved Drug")
            })

df = pd.DataFrame(records)
print(f"\nRebuilt raw dataset with {len(df)} rows.")

# 1. Deduplicate by chembl_id in case of ChEMBL API weirdness
df = df.drop_duplicates(subset="chembl_id")

# 2. Deduplicate by drug_name
df = df.drop_duplicates(subset="drug_name")
print(f"After drug_name dedup: {len(df)}")

# 3. Deduplicate by 25-character InChIKey (Skeleton + Stereochem)
def get_inchi_25(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    inchi = Chem.MolToInchiKey(mol)
    if not inchi: return None
    return inchi[:25]

df["inchikey_25"] = df["smiles"].apply(get_inchi_25)

before = len(df)
df = df.drop_duplicates(subset="inchikey_25", keep="first")
df = df.drop(columns=["inchikey_25"])
after = len(df)

print(f"After 25-char InChIKey dedup: {after} (Dropped {before - after} true structural duplicates/salts)")
print(f"This is {156 - (before - after)} FEWER duplicates than the 14-char dedup dropped.")

df.to_csv("full_dataset.csv", index=False)
print("Saved to full_dataset.csv")

import build_db
build_db.build_database()
