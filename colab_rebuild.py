# Run this script in a Google Colab notebook.
# First, install RDKit in the notebook by running: !pip install rdkit pandas
# Upload the 'full_dataset.csv' to the Colab environment before running.

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser, Uncharger
from rdkit.Chem import DataStructs
import json

def run_rebuild():
    print("Loading raw compound pool...")
    df = pd.read_csv("full_dataset.csv")
    print(f"Raw pull count: {len(df)}")
    
    lfc = LargestFragmentChooser()
    uncharger = Uncharger()
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    
    records = []
    fp_objects = []  # Store native RDKit fingerprint objects for fast hubness calc
    failed_rows = [] # Track explicit failures
    
    for i, row in df.iterrows():
        smiles = row['smiles']
        if pd.isna(smiles) or not smiles:
            continue
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
            
        # 1. Desalt using LargestFragmentChooser and Uncharger
        try:
            mol = lfc.choose(mol)
            mol = uncharger.uncharge(mol)
        except Exception as e:
            failed_rows.append({"drug_name": row['drug_name'], "error": str(e), "step": "cleaning"})
            continue
            
        # 2. Resanitize
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            failed_rows.append({"drug_name": row['drug_name'], "error": str(e), "step": "sanitize"})
            continue
            
        # 3. Generate Canonical SMILES and InChIKey
        canonical = Chem.MolToSmiles(mol)
        if not canonical: continue
            
        inchikey = Chem.MolToInchiKey(mol)
        
        # 4. Generate Chirality-Aware Morgan Fingerprint
        fp = morgan_gen.GetFingerprint(mol)
        fp_bitstring = fp.ToBitString()
        
        records.append({
            "drug_name": row['drug_name'],
            "chembl_id": row['chembl_id'],
            "original_smiles": smiles,
            "desalted_smiles": canonical,
            "inchikey": inchikey,
            "atc_category": row['atc_category'],
            "indication": row['indication'],
            "morgan_fp": fp_bitstring
        })
        fp_objects.append(fp)
        
    df_processed = pd.DataFrame(records)
    
    # Associate native FP objects securely with rows so dedup drops the correct ones
    df_processed['native_fp'] = fp_objects
    
    # 5. Deduplicate by InChIKey (keep first)
    df_dedup = df_processed.drop_duplicates(subset="inchikey", keep="first").copy()
    print(f"Post-desalting and deduplication count: {len(df_dedup)}")
    
    # Check for residual multi-component structures (Vardenafil-style tautomer issues)
    multi_comp = df_dedup['desalted_smiles'].str.contains('.', regex=False).sum()
    if multi_comp > 0:
        print(f"WARNING: Found {multi_comp} structures still containing multiple components ('.') after LargestFragmentChooser.")
    else:
        print("SUCCESS: No multi-component structures remaining.")
    
    # Validation Checks
    benchmarks = ['SILDENAFIL', 'MINOXIDIL', 'METFORMIN', 'THALIDOMIDE', 'DULOXETINE', 'GABAPENTIN', 'BUPROPION', 'RALOXIFENE']
    missing = [d for d in benchmarks if not df_dedup['drug_name'].str.contains(d).any()]
    if not missing:
        print("SUCCESS: All 8 benchmark drugs survived the dedup pass.")
    else:
        print(f"WARNING: Missing benchmarks: {missing}")
        
    # Chirality Check
    ofloxacin = df_dedup[df_dedup['drug_name'].str.contains('OFLOXACIN', na=False)]
    levofloxacin = df_dedup[df_dedup['drug_name'].str.contains('LEVOFLOXACIN', na=False)]
    
    if not ofloxacin.empty and not levofloxacin.empty:
        oflox_fp = ofloxacin.iloc[0]['morgan_fp']
        levoflox_fp = levofloxacin.iloc[0]['morgan_fp']
        
        if oflox_fp != levoflox_fp:
            print("CHIRALITY CHECK PASSED: Ofloxacin and Levofloxacin have distinct Morgan fingerprints.")
        else:
            print("CHIRALITY CHECK FAILED: Ofloxacin and Levofloxacin fingerprints are identical. (Is includeChirality=True working?)")
    else:
        print("CHIRALITY CHECK INCONCLUSIVE: Missing one or both of the drugs.")
        
    print("Computing hubness scores (average Tanimoto similarity to all other drugs)...")
    
    # Extract the deduplicated native FP objects
    final_fp_list = df_dedup['native_fp'].tolist()
    
    hubness_scores = []
    for i in range(len(final_fp_list)):
        sims = DataStructs.BulkTanimotoSimilarity(final_fp_list[i], final_fp_list)
        # Average similarity (excluding self-similarity of 1.0)
        avg_sim = (sum(sims) - 1.0) / (len(final_fp_list) - 1) if len(final_fp_list) > 1 else 0
        hubness_scores.append(avg_sim)
        
    df_dedup['hubness_score'] = hubness_scores
    df_dedup['hubness_percentile'] = df_dedup['hubness_score'].rank(pct=True) * 100
    
    # Drop the temporary native_fp column before export
    df_dedup = df_dedup.drop(columns=['native_fp'])
    
    # Export the final dataset
    if failed_rows:
        print(f"\nWARNING: {len(failed_rows)} rows explicitly failed the cleaning/sanitization pipeline:")
        for fail in failed_rows[:10]:
            print(f"  - {fail['drug_name']} ({fail['step']}): {fail['error']}")
        if len(failed_rows) > 10:
            print("  ... and more.")
    else:
        print("\nSUCCESS: 0 rows explicitly failed the pipeline.")
        
    df_dedup.to_csv("colab_output_dataset.csv", index=False)
    print("Exported to colab_output_dataset.csv")

if __name__ == "__main__":
    run_rebuild()
