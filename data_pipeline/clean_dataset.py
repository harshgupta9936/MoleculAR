import pandas as pd
import json
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def clean_existing_dataset():
    print("Loading full_dataset.csv...")
    df = pd.read_csv("full_dataset.csv")
    
    largest_fragment = rdMolStandardize.LargestFragmentChooser()
    
    def process_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        try:
            mol = largest_fragment.choose(mol)
        except:
            pass
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)
    
    print(f"Cleaning {len(df)} SMILES strings...")
    df["clean_smiles"] = df["smiles"].apply(process_smiles)
    
    df = df.dropna(subset=["clean_smiles"])
    df["smiles"] = df["clean_smiles"]
    df = df.drop(columns=["clean_smiles"])
    
    before_dedup = len(df)
    df = df.drop_duplicates(subset="smiles", keep="first")
    after_dedup = len(df)
    
    print(f"Deduplication removed {before_dedup - after_dedup} identical active fragments (salts/hydrates).")
    print(f"New dataset size: {len(df)}")
    
    # Save the cleaned CSV
    df.to_csv("full_dataset.csv", index=False)
    
    # Filter the target json
    valid_ids = set(df["chembl_id"])
    with open("drug_targets.json", "r") as f:
        targets = json.load(f)
    
    cleaned_targets = {k: v for k, v in targets.items() if k in valid_ids}
    with open("drug_targets.json", "w") as f:
        json.dump(cleaned_targets, f, indent=2)
        
    print("Cleaned up dataset and aligned JSONs successfully!")

if __name__ == "__main__":
    clean_existing_dataset()
