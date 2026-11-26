import pandas as pd
import sqlite3
import re
from rdkit import DataStructs

def jaccard_similarity(set1, set2):
    set1, set2 = set(set1), set(set2)
    if not set1 or not set2:
        return 0.0
    return len(set1 & set2) / len(set1 | set2)

def is_valid_drug_match(iname, dname):
    return bool(re.search(rf'\b{re.escape(dname)}\b', str(iname), re.IGNORECASE))

def load_drug_features():
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    cursor.execute("SELECT chembl_id, generic_base_name, fingerprint_blob, hubness_percentile FROM drugs WHERE generic_base_name IS NOT NULL")
    
    drugs = {}
    
    for cid, name, fp_blob, hub_pct in cursor.fetchall():
        dname = str(name).strip().upper()
        if not dname: continue
        fp = None
        if fp_blob:
            try:
                fp = DataStructs.cDataStructs.CreateFromBitString(str(fp_blob))
            except:
                pass
            
        drugs[dname] = {
            'fingerprint': fp,
            'hubness_percentile': hub_pct if hub_pct is not None else 50.0,
            'targets': set(),
            'side_effects': set()
        }
        
    cursor.execute("SELECT d.generic_base_name, t.uniprot_id, t.target_name FROM targets t JOIN drugs d ON t.chembl_id = d.chembl_id")
    for name, uid, tname in cursor.fetchall():
        dname = str(name).strip().upper()
        if dname in drugs: 
            val = uid if uid else tname
            if val:
                drugs[dname]['targets'].add(val.lower())
        
    cursor.execute("SELECT d.generic_base_name, s.side_effect_name FROM side_effects s JOIN drugs d ON s.chembl_id = d.chembl_id")
    for name, se in cursor.fetchall():
        dname = str(name).strip().upper()
        if dname in drugs: 
            if se:
                drugs[dname]['side_effects'].add(se.lower())
        
    conn.close()
    return drugs

def get_disease_references(drug_names):
    print("Parsing AACT for reference drugs...")
    cond_df = pd.read_csv('conditions.txt', sep='|', usecols=['nct_id', 'downcase_name'], on_bad_lines='skip', low_memory=False)
    int_df = pd.read_csv('interventions.txt', sep='|', usecols=['nct_id', 'intervention_type', 'name'], on_bad_lines='skip', low_memory=False)
    int_df = int_df[int_df['intervention_type'].astype(str).str.upper().isin(['DRUG', 'BIOLOGICAL', 'DIETARY SUPPLEMENT'])]
    
    unique_interventions = int_df['name'].dropna().unique()
    intervention_mapping = {}
    
    drug_names_lower = {d.lower(): d for d in drug_names}
    sorted_drugs = sorted(drug_names_lower.keys(), key=len, reverse=True)
    pattern = re.compile(r'\b(?:' + '|'.join(map(re.escape, sorted_drugs)) + r')\b', re.IGNORECASE)
    
    for iname in unique_interventions:
        found = pattern.findall(str(iname))
        if found:
            matched = list({drug_names_lower[f.lower()] for f in found if f.lower() in drug_names_lower})
            if matched:
                intervention_mapping[iname] = matched
            
    nct_to_drugs = {}
    for row in int_df.itertuples(index=False):
        nid = row.nct_id
        iname = row.name
        matched_drugs = intervention_mapping.get(iname, [])
        for dname in matched_drugs:
            if nid not in nct_to_drugs: nct_to_drugs[nid] = set()
            nct_to_drugs[nid].add(dname)
            
    nct_to_cond = {}
    for row in cond_df.itertuples(index=False):
        nid = row.nct_id
        cname = str(row.downcase_name).lower()
        if cname == 'nan': continue
        if nid not in nct_to_cond: nct_to_cond[nid] = set()
        nct_to_cond[nid].add(cname)
        
    cond_to_drugs = {}
    for nid, dnames in nct_to_drugs.items():
        if nid in nct_to_cond:
            for cname in nct_to_cond[nid]:
                if cname not in cond_to_drugs: cond_to_drugs[cname] = set()
                cond_to_drugs[cname].update(dnames)
                
    return cond_to_drugs

def main():
    print("Loading drug features...")
    drugs = load_drug_features()
    
    cond_to_drugs = get_disease_references(list(drugs.keys()))
    
    df = pd.read_csv('training_labels_v2.csv')
    print(f"Loaded {len(df)} pairs. Computing features...")
    
    records = []
    
    for idx, row in df.iterrows():
        drug = str(row['drug_name']).upper()
        disease = str(row['disease']).lower()
        label = row['label']
        
        q_feats = drugs.get(drug)
        if not q_feats:
            print(f"Missing features for {drug}")
            continue
            
        ref_drugs = cond_to_drugs.get(disease, set()) - {drug}
        
        r_fps = []
        agg_targets = set()
        agg_ses = set()
        
        for rd in ref_drugs:
            r_feat = drugs.get(rd)
            if not r_feat: continue
            if r_feat['fingerprint']:
                r_fps.append(r_feat['fingerprint'])
            agg_targets.update(r_feat['targets'])
            agg_ses.update(r_feat['side_effects'])
            
        max_struct = 0.0
        if q_feats['fingerprint'] and r_fps:
            sims = DataStructs.BulkTanimotoSimilarity(q_feats['fingerprint'], r_fps)
            max_struct = max(sims)
            
        target_sim = jaccard_similarity(q_feats['targets'], agg_targets)
        se_sim = jaccard_similarity(q_feats['side_effects'], agg_ses)
        hubness = q_feats['hubness_percentile']
        
        records.append({
            'drug_name': drug,
            'disease': disease,
            'label': label,
            'struct_sim': max_struct,
            'target_sim': target_sim,
            'side_effect_sim': se_sim,
            'hubness_percentile': hubness
        })
        
    out_df = pd.DataFrame(records)
    out_df.to_csv('feature_matrix.csv', index=False)
    print(f"Generated feature matrix with {len(out_df)} rows. Saved to feature_matrix.csv.")
    
if __name__ == "__main__":
    main()
