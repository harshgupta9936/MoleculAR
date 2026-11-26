"""
Recomputes hubness scores and percentile ranks for all valid (non-trivial) drugs in drugrepo.db.

Methodology Note & Percentile Scope:
- In accordance with methodology notes (Task 2.4), trivial compounds (gases, radioisotopes,
  simple inorganic salts, and tiny fragments) are excluded from the hubness distribution
  so they do not artificially distort percentile ranks for genuine drug molecules.
- Non-trivial compounds (2,202) have their average pairwise Tanimoto similarity computed
  strictly within the non-trivial drug pool, and their percentiles ranked from (1 / N_valid) * 100
  up to 100.0%.
- Trivial compounds (65) have their hubness_score and hubness_percentile set to NULL.
- The pipeline (Mode A and Mode B in pipeline.py) safely handles NULL hubness without TypeError.
"""

import sqlite3
import pandas as pd
from rdkit.Chem import DataStructs

def recompute_hubness(db_path='drugrepo.db'):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # 1. Fetch valid (non-trivial) fingerprints
    cursor.execute("""
        SELECT chembl_id, fingerprint_blob 
        FROM drugs 
        WHERE trivial_category IS NULL AND fingerprint_blob IS NOT NULL
    """)
    rows = cursor.fetchall()
    
    print(f"Recomputing hubness for {len(rows)} non-trivial drugs...")
    
    cids = []
    fps = []
    for cid, blob in rows:
        cids.append(cid)
        fps.append(DataStructs.CreateFromBitString(blob))
        
    hubness_scores = []
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
        avg_sim = (sum(sims) - 1.0) / (len(fps) - 1) if len(fps) > 1 else 0.0
        hubness_scores.append(avg_sim)
        
    # Calculate percentiles strictly within the non-trivial drug pool
    df = pd.DataFrame({"chembl_id": cids, "hubness_score": hubness_scores})
    df['hubness_percentile'] = df['hubness_score'].rank(pct=True) * 100.0
    
    # Update DB for non-trivial compounds
    updates = []
    for _, row in df.iterrows():
        updates.append((float(row['hubness_score']), float(row['hubness_percentile']), row['chembl_id']))
        
    cursor.executemany("UPDATE drugs SET hubness_score=?, hubness_percentile=? WHERE chembl_id=?", updates)
    conn.commit()
    print(f"Updated hubness scores and percentiles for {len(updates)} non-trivial drugs.")
    
    # Explicitly clear hubness for trivial compounds per methodology notes
    cursor.execute("UPDATE drugs SET hubness_score=NULL, hubness_percentile=NULL WHERE trivial_category IS NOT NULL")
    cleared_count = cursor.rowcount
    conn.commit()
    print(f"Cleared hubness to NULL for {cleared_count} trivial compounds.")
    
    conn.close()

if __name__ == "__main__":
    recompute_hubness()
