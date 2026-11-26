import sqlite3
from rdkit import Chem
from rdkit.Chem import Descriptors

def main():
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    
    try:
        cursor.execute("ALTER TABLE drugs ADD COLUMN trivial_category TEXT")
    except:
        pass
        
    cursor.execute("SELECT chembl_id, drug_name, smiles, generic_base_name FROM drugs")
    all_drugs = cursor.fetchall()
    
    trivial_candidates = {}
    
    for row in all_drugs:
        cid, name, smiles, base = row
        if not smiles: continue
        
        is_trivial = False
        cat = ""
        
        # Determine category
        name_upper = name.upper()
        if any(x in name_upper for x in [" 89", " 127", " 133", " 82", " 18", " 123", " 131", " 81M", " 13", " 32", " 223", " 129"]):
            # Fast heuristic for Radioisotopes in the names
            is_trivial = True
            cat = "Radioisotope"
        elif any(x in name_upper for x in ["XENON", "KRYPTON", "HELIUM", "NITROGEN", "OXYGEN", "NITROUS OXIDE", "NITRIC OXIDE"]):
            is_trivial = True
            cat = "Medical Gas"
        elif 'c' not in smiles.lower():
            is_trivial = True
            if "AMMONIA" in name_upper:
                cat = "Medical Gas"
            else:
                cat = "Simple Salt/Inorganic"
        else:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    hw = Descriptors.HeavyAtomCount(mol)
                    mw = Descriptors.MolWt(mol)
                    if mw < 60 or hw < 4:
                        is_trivial = True
                        cat = "Tiny Organic"
            except:
                pass
                
        if is_trivial:
            if not cat:
                cat = "Trivial/Inorganic"
            trivial_candidates[cid] = (name, cat)
            cursor.execute("UPDATE drugs SET trivial_category = ? WHERE chembl_id = ?", (cat, cid))
    
    conn.commit()
    
    # Assert no overlap with benchmark/hubs
    benchmarks = {"SILDENAFIL", "MINOXIDIL", "METFORMIN", "THALIDOMIDE", "DULOXETINE", "GABAPENTIN", "BUPROPION", "RALOXIFENE"}
    hubs = {"CLOZAPINE", "OLANZAPINE", "QUETIAPINE", "NEFAZODONE", "FENOFIBRATE", "DICLOFENAC"}
    
    overlaps = []
    for cid, (name, cat) in trivial_candidates.items():
        if name in benchmarks or name in hubs:
            overlaps.append(name)
            
    print(f"Tagged {len(trivial_candidates)} trivial candidates.")
    print(f"Overlap with locked benchmark/hub drugs: {len(overlaps)}")
    if overlaps:
        print(overlaps)
    conn.close()

if __name__ == "__main__":
    main()
