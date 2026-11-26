import sqlite3
import gzip
import urllib.request
import json
import random
from rdkit import Chem
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser, Uncharger

def standardize_smiles(smiles):
    if not smiles or not isinstance(smiles, str):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
            
        # 1. Largest Fragment Chooser
        lfc = LargestFragmentChooser()
        mol = lfc.choose(mol)
        
        # 2. Uncharger
        uc = Uncharger()
        mol = uc.uncharge(mol)
        
        return Chem.MolToInchiKey(mol)
    except Exception as e:
        print("Error standardizing:", e)
        return None

def fetch_pubchem_smiles(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES,CanonicalSMILES/JSON"
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=5) as response:
            data = json.loads(response.read().decode('utf-8'))
            props = data['PropertyTable']['Properties'][0]
            if 'IsomericSMILES' in props:
                return props['IsomericSMILES']
            elif 'CanonicalSMILES' in props:
                return props['CanonicalSMILES']
            else:
                for k, v in props.items():
                    if 'SMILES' in k:
                        return v
            return None
    except Exception as e:
        return None

def main():
    print("Building Crosswalk...")
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    
    # Create crosswalk table
    cursor.execute("DROP TABLE IF EXISTS crosswalk")
    cursor.execute("""
        CREATE TABLE crosswalk (
            chembl_id TEXT PRIMARY KEY,
            drugbank_id TEXT,
            pubchem_cid TEXT,
            FOREIGN KEY(chembl_id) REFERENCES drugs(chembl_id)
        )
    """)
    
    # Get all 2,267 chembl_ids
    cursor.execute("SELECT chembl_id, smiles FROM drugs")
    drugs_data = cursor.fetchall()
    target_cids = {row[0]: {'smiles': row[1], 'drugbank_id': None, 'pubchem_cid': None} for row in drugs_data}
    
    # Parse src1src2.txt.gz (ChEMBL -> DrugBank)
    print("Parsing ChEMBL to DrugBank...")
    with gzip.open('unichem_data/src1src2.txt.gz', 'rt') as f:
        f.readline() # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cid, dbid = parts
                if cid in target_cids:
                    target_cids[cid]['drugbank_id'] = dbid
                    
    # Parse src1src22.txt.gz (ChEMBL -> PubChem)
    print("Parsing ChEMBL to PubChem...")
    with gzip.open('unichem_data/src1src22.txt.gz', 'rt') as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                cid, pcid = parts
                if cid in target_cids:
                    target_cids[cid]['pubchem_cid'] = pcid
                    
    # Insert into DB
    for cid, ids in target_cids.items():
        cursor.execute("INSERT INTO crosswalk (chembl_id, drugbank_id, pubchem_cid) VALUES (?, ?, ?)",
                      (cid, ids['drugbank_id'], ids['pubchem_cid']))
    conn.commit()
    
    # Stats
    cursor.execute("SELECT COUNT(*) FROM crosswalk")
    total = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM crosswalk WHERE drugbank_id IS NOT NULL")
    db_count = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM crosswalk WHERE pubchem_cid IS NOT NULL")
    pc_count = cursor.fetchone()[0]
    
    print(f"\nCrosswalk Table Row Count: {total}")
    print(f"Mapped to DrugBank: {db_count} ({(db_count/total)*100:.1f}%)")
    print(f"Mapped to PubChem: {pc_count} ({(pc_count/total)*100:.1f}%)")
    
    # Validation
    print("\nValidating 50 Random Drugs...")
    cursor.execute("SELECT chembl_id, drugbank_id, pubchem_cid, smiles FROM drugs JOIN crosswalk USING (chembl_id) WHERE pubchem_cid IS NOT NULL AND drugbank_id IS NOT NULL")
    validatable = cursor.fetchall()
    
    sample = random.sample(validatable, min(50, len(validatable)))
    mismatches = 0
    success = 0
    
    for row in sample:
        cid, dbid, pcid, c_smiles = row
        p_smiles = fetch_pubchem_smiles(pcid)
        if not p_smiles:
            print(f"Failed to fetch PubChem SMILES for {pcid}")
            continue
            
        c_inchi = standardize_smiles(c_smiles)
        p_inchi = standardize_smiles(p_smiles)
        
        if c_inchi and p_inchi:
            if c_inchi == p_inchi:
                success += 1
            else:
                mismatches += 1
                print(f"MISMATCH: {cid} (DB: {dbid}, PC: {pcid})")
                print(f"  ChEMBL InChIKey:  {c_inchi}")
                print(f"  PubChem InChIKey: {p_inchi}")
        else:
            print(f"Standardization failed for {cid}")
            
    print(f"\nValidation complete.")
    print(f"Tested: 50, Matched: {success}, Mismatches: {mismatches}")
    
    conn.close()

if __name__ == "__main__":
    main()
