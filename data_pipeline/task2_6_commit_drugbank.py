import urllib.request
import pandas as pd
import io
import sqlite3
import requests
import urllib.parse
import concurrent.futures
from rdkit import Chem
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser, Uncharger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import Descriptors
import time

def standardize_smiles_and_mol(mol):
    try:
        if not mol: return None, None
        lfc = LargestFragmentChooser()
        mol = lfc.choose(mol)
        uc = Uncharger()
        mol = uc.uncharge(mol)
        return Chem.MolToInchiKey(mol), mol
    except:
        return None, None

def check_openfda(base_name):
    url = f'https://api.fda.gov/drug/label.json?search=openfda.generic_name:"{urllib.parse.quote(base_name)}"&limit=1'
    try:
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            return True
    except:
        pass
    return False

def check_dailymed(base_name):
    url = f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?drug_name={urllib.parse.quote(base_name)}"
    try:
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            data = r.json()
            if data.get("data") and len(data["data"]) > 0:
                return True
    except:
        pass
    return False

def audit_market_status(item):
    cid, name, groups = item
    time.sleep(0.05)
    
    vet = 'vet_approved' in groups
    withdrawn = 'withdrawn' in groups
    
    if vet:
        return cid, 'Veterinary'
    if withdrawn:
        return cid, 'Withdrawn'
        
    if check_openfda(name) or check_dailymed(name):
        return cid, 'Active (US)'
        
    return cid, 'International/Non-US'

def main():
    print("Downloading DrugBank Open TSV...")
    url = 'https://raw.githubusercontent.com/dhimmel/drugbank/master/data/drugbank.tsv'
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as r:
        df = pd.read_csv(io.StringIO(r.read().decode('utf-8')), sep='\t')
        
    approved = df[df['groups'].str.contains('approved', na=False) & (df['type'] == 'small molecule')]
    
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    
    cursor.execute("SELECT chembl_id, smiles, fingerprint_blob FROM drugs")
    existing = cursor.fetchall()
    existing_inchis = set()
    for cid, smiles, blob in existing:
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        ikey, _ = standardize_smiles_and_mol(mol)
        if ikey:
            existing_inchis.add(ikey)
            
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    
    new_drugs = []
    
    for i, row in approved.iterrows():
        inchi = row.get('inchi')
        mol = None
        if pd.notna(inchi):
            mol = Chem.MolFromInchi(inchi)
            
        ikey, st_mol = standardize_smiles_and_mol(mol)
        if not ikey: continue
            
        if ikey not in existing_inchis:
            existing_inchis.add(ikey)
            db_id = row['drugbank_id']
            name = str(row['name']).upper()
            smiles = Chem.MolToSmiles(st_mol)
            
            # Trivial check
            is_trivial = False
            cat = None
            if 'c' not in smiles.lower():
                is_trivial = True
                cat = "Simple Salt/Inorganic"
            else:
                hw = Descriptors.HeavyAtomCount(st_mol)
                mw = Descriptors.MolWt(st_mol)
                if mw < 60 or hw < 4:
                    is_trivial = True
                    cat = "Tiny Organic"
                    
            fp = morgan_gen.GetFingerprint(st_mol)
            fp_bitstring = fp.ToBitString()
            
            new_drugs.append({
                'chembl_id': db_id,
                'drug_name': name,
                'smiles': smiles,
                'generic_base_name': name,
                'trivial_category': cat,
                'groups': str(row['groups']),
                'fingerprint_blob': fp_bitstring
            })
            
    print(f"Inserting {len(new_drugs)} new DrugBank drugs...")
    inserts = []
    for d in new_drugs:
        inserts.append((d['chembl_id'], d['drug_name'], d['smiles'], 'Active (US)', d['generic_base_name'], 'DrugBank', d['fingerprint_blob'], d['trivial_category']))
        
    cursor.executemany("""
        INSERT INTO drugs (chembl_id, drug_name, smiles, market_status, generic_base_name, source, fingerprint_blob, trivial_category)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """, inserts)
    conn.commit()
    
    print("Running market status audit via openFDA/DailyMed...")
    audit_tasks = [(d['chembl_id'], d['generic_base_name'], d['groups']) for d in new_drugs]
    results = {'Active (US)': 0, 'Veterinary': 0, 'Withdrawn': 0, 'International/Non-US': 0}
    
    updates = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        for res in executor.map(audit_market_status, audit_tasks):
            cid, market_status = res
            updates.append((market_status, cid))
            results[market_status] += 1
            
    cursor.executemany("UPDATE drugs SET market_status = ? WHERE chembl_id = ?", updates)
    conn.commit()
    
    for k, v in results.items():
        print(f"  {k}: {v}")
        
    # Trivial count
    trivial = sum(1 for d in new_drugs if d['trivial_category'])
    print(f"Trivial/Inorganic compounds flagged in new batch: {trivial}")
    
    conn.close()

if __name__ == "__main__":
    main()
