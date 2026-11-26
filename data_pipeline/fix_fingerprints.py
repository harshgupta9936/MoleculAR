import sqlite3
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

def main():
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    
    cursor.execute("SELECT chembl_id, smiles FROM drugs WHERE fingerprint_blob IS NULL AND trivial_category IS NULL")
    rows = cursor.fetchall()
    print(f"Fixing {len(rows)} drugs with missing fingerprints...")
    
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    
    updates = []
    for cid, smiles in rows:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = morgan_gen.GetFingerprint(mol)
            fp_bitstring = fp.ToBitString()
            updates.append((fp_bitstring, cid))
            
    if updates:
        cursor.executemany("UPDATE drugs SET fingerprint_blob = ? WHERE chembl_id = ?", updates)
        conn.commit()
        print(f"Updated {len(updates)} fingerprints.")
    
    conn.close()

if __name__ == "__main__":
    main()
