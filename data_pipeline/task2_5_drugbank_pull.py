import urllib.request
import pandas as pd
import io
import sqlite3
from rdkit import Chem
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser, Uncharger

def standardize_smiles(mol):
    try:
        if not mol: return None
        lfc = LargestFragmentChooser()
        mol = lfc.choose(mol)
        uc = Uncharger()
        mol = uc.uncharge(mol)
        return Chem.MolToInchiKey(mol)
    except:
        return None

def main():
    print("Downloading DrugBank Open TSV...")
    url = 'https://raw.githubusercontent.com/dhimmel/drugbank/master/data/drugbank.tsv'
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as r:
        df = pd.read_csv(io.StringIO(r.read().decode('utf-8')), sep='\t')
        
    approved = df[df['groups'].str.contains('approved', na=False) & (df['type'] == 'small molecule')]
    print(f"Found {len(approved)} approved small molecules in DrugBank Open Data.")
    
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    
    cursor.execute("SELECT chembl_id, smiles, fingerprint_blob FROM drugs")
    existing = cursor.fetchall()
    existing_inchis = set()
    for cid, smiles, blob in existing:
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        ikey = standardize_smiles(mol)
        if ikey:
            existing_inchis.add(ikey)
            
    print(f"Loaded {len(existing_inchis)} valid existing InChIKeys from our 2290 pool.")
    
    new_drugs = 0
    collisions = 0
    
    for i, row in approved.iterrows():
        # Try to parse from InChI first since it's most robust
        inchi = row.get('inchi')
        mol = None
        if pd.notna(inchi):
            mol = Chem.MolFromInchi(inchi)
            
        ikey = standardize_smiles(mol)
        if not ikey:
            continue
            
        if ikey in existing_inchis:
            collisions += 1
        else:
            new_drugs += 1
            existing_inchis.add(ikey) # Prevent internal duplicates
            
    print(f"DrugBank -> Existing Collisions: {collisions}")
    print(f"DrugBank -> Genuinely New Drugs: {new_drugs}")
    
    conn.close()

if __name__ == "__main__":
    main()
