import sqlite3
import urllib.request
import urllib.parse
import json
import time
from rdkit import Chem
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser, Uncharger

def fetch_pubchem_bulk(cids):
    # Fetch in batches of 500
    results = {}
    batch_size = 500
    for i in range(0, len(cids), batch_size):
        batch = cids[i:i+batch_size]
        data = urllib.parse.urlencode({'cid': ','.join(batch)}).encode('utf-8')
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/IsomericSMILES,CanonicalSMILES/JSON"
        req = urllib.request.Request(url, data=data)
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                resp_json = json.loads(response.read().decode('utf-8'))
                for props in resp_json['PropertyTable']['Properties']:
                    cid = str(props['CID'])
                    smiles = props.get('IsomericSMILES') or props.get('CanonicalSMILES') or props.get('SMILES') or props.get('ConnectivitySMILES')
                    results[cid] = smiles
        except Exception as e:
            print(f"Batch {i} failed: {e}")
        time.sleep(1) # Be polite
    return results

def standardize_smiles(smiles):
    if not smiles or not isinstance(smiles, str):
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        lfc = LargestFragmentChooser()
        mol = lfc.choose(mol)
        uc = Uncharger()
        mol = uc.uncharge(mol)
        return Chem.MolToInchiKey(mol)
    except:
        return None

def main():
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    
    # 1. Add source column if missing
    try:
        cursor.execute("ALTER TABLE drugs ADD COLUMN source VARCHAR(100) DEFAULT 'ChEMBL'")
        conn.commit()
    except:
        pass
        
    # Generate baseline InChIKeys for ChEMBL (in memory)
    cursor.execute("SELECT chembl_id, smiles, source FROM drugs")
    chembl_drugs = cursor.fetchall()
    existing_inchis = {}
    for cid, smiles, src in chembl_drugs:
        ikey = standardize_smiles(smiles)
        if ikey:
            existing_inchis[ikey] = cid
            
    print(f"Baseline ChEMBL drugs with valid InChIKeys: {len(existing_inchis)}")
            
    # 2. Get mapped PubChem CIDs
    cursor.execute("SELECT chembl_id, pubchem_cid FROM crosswalk WHERE pubchem_cid IS NOT NULL")
    mapped_pc = cursor.fetchall()
    cids_to_fetch = [row[1] for row in mapped_pc]
    cid_to_chembl = {row[1]: row[0] for row in mapped_pc}
    
    print(f"Fetching {len(cids_to_fetch)} PubChem structures via bulk POST...")
    pc_smiles_dict = fetch_pubchem_bulk(cids_to_fetch)
    print(f"Fetched {len(pc_smiles_dict)} SMILES from PubChem.")
    
    # 3. Process PubChem
    pc_collisions_existing = 0
    pc_collisions_internal = 0
    pc_new = 0
    pc_processed_inchis = set()
    
    mismatched_chembl_mappings = []
    
    for pcid, smiles in pc_smiles_dict.items():
        ikey = standardize_smiles(smiles)
        if not ikey: continue
        
        expected_chembl = cid_to_chembl.get(pcid)
        
        if ikey in existing_inchis:
            pc_collisions_existing += 1
            matched_chembl = existing_inchis[ikey]
            if matched_chembl != expected_chembl:
                mismatched_chembl_mappings.append((pcid, expected_chembl, matched_chembl))
            # Tag existing with PubChem
            cursor.execute("UPDATE drugs SET source = source || ', PubChem' WHERE chembl_id = ? AND source NOT LIKE '%PubChem%'", (matched_chembl,))
        elif ikey in pc_processed_inchis:
            pc_collisions_internal += 1
        else:
            pc_new += 1
            pc_processed_inchis.add(ikey)
            
    conn.commit()
    print(f"PubChem -> Existing Collisions: {pc_collisions_existing}")
    print(f"PubChem -> Internal Collisions: {pc_collisions_internal}")
    print(f"PubChem -> New Unique Drugs: {pc_new}")
    if mismatched_chembl_mappings:
        print(f"WARNING: {len(mismatched_chembl_mappings)} PubChem structures matched a DIFFERENT ChEMBL drug than UniChem predicted based on strict InChIKey!")
        
    conn.close()

if __name__ == "__main__":
    main()
