import pandas as pd
import sqlite3
import random
import re
import json

def is_circular(cond, ind):
    cond = cond.lower()
    ind = ind.lower()
    if cond in ind: return True
    
    synonyms = [
        ('cancer', 'neoplasm'), ('cancer', 'carcinoma'), ('cancer', 'tumor'), ('cancer', 'malignan'), ('cancer', 'sarcoma'),
        ('heart', 'cardiac'), ('heart', 'myocardial'),
        ('kidney', 'renal'), ('liver', 'hepatic'), ('lung', 'pulmonary'),
        ('brain', 'cerebral'), ('stomach', 'gastric'), ('blood', 'hematologic')
    ]
    
    for s1, s2 in synonyms:
        if s1 in cond and s2 in ind: return True
        if s2 in cond and s1 in ind: return True
        
    words1 = set(re.findall(r'\b\w{5,}\b', cond))
    words2 = set(re.findall(r'\b\w{5,}\b', ind))
    if words1 & words2: return True
    
    return False

def main():
    print("Loading drugs...")
    conn = sqlite3.connect('drugrepo.db')
    c = conn.cursor()
    c.execute("SELECT chembl_id, generic_base_name, indication FROM drugs WHERE indication != 'Other' AND indication IS NOT NULL AND indication != '' AND trivial_category IS NULL")
    
    drugs = {}
    for cid, name, ind in c.fetchall():
        dname = str(name).strip().upper()
        drugs[dname] = str(ind).lower()
    conn.close()
    
    print(f"Loaded {len(drugs)} active drugs with real indications.")
    
    print("Filtering studies...")
    studies = pd.read_csv('studies.txt', sep='|', usecols=['nct_id', 'study_type', 'phase'], on_bad_lines='skip', low_memory=False)
    studies = studies[(studies['study_type'].astype(str).str.upper() == 'INTERVENTIONAL') & 
                      (studies['phase'].astype(str).str.upper().isin(['PHASE2', 'PHASE3', 'PHASE4']))]
    valid_ncts = set(studies['nct_id'].unique())
    print(f"Found {len(valid_ncts)} valid interventional phase 2-4 studies.")
    
    print("Filtering designs...")
    designs = pd.read_csv('designs.txt', sep='|', usecols=['nct_id', 'primary_purpose'], on_bad_lines='skip', low_memory=False)
    designs = designs[designs['primary_purpose'].astype(str).str.upper() == 'TREATMENT']
    treatment_ncts = set(designs['nct_id'].unique())
    
    valid_ncts = valid_ncts & treatment_ncts
    print(f"Found {len(valid_ncts)} studies primarily for TREATMENT.")
    
    print("Mapping conditions...")
    cond_df = pd.read_csv('conditions.txt', sep='|', usecols=['nct_id', 'downcase_name'], on_bad_lines='skip', low_memory=False)
    cond_df = cond_df[cond_df['nct_id'].isin(valid_ncts)]
    
    nct_to_cond = {}
    for row in cond_df.itertuples(index=False):
        nid = row.nct_id
        cname = str(row.downcase_name).lower()
        if cname == 'nan': continue
        if nid not in nct_to_cond: nct_to_cond[nid] = set()
        nct_to_cond[nid].add(cname)
        
    print("Mapping interventions...")
    int_df = pd.read_csv('interventions.txt', sep='|', usecols=['nct_id', 'intervention_type', 'name'], on_bad_lines='skip', low_memory=False)
    int_df = int_df[(int_df['nct_id'].isin(valid_ncts)) & (int_df['intervention_type'].astype(str).str.upper().isin(['DRUG', 'BIOLOGICAL', 'DIETARY SUPPLEMENT']))]
    
    drug_names = list(drugs.keys())
    drug_names.sort(key=len, reverse=True)
    int_df['name_upper'] = int_df['name'].astype(str).str.upper()
    unique_interventions = int_df['name_upper'].unique()
    
    intervention_mapping = {}
    for iname in unique_interventions:
        if not isinstance(iname, str): continue
        matched = []
        for dname in drug_names:
            if len(dname) < 5: continue
            if dname in iname:
                matched.append(dname)
        if matched:
            intervention_mapping[iname] = matched
            
    drug_to_conds = {dname: set() for dname in drug_names}
    for row in int_df.itertuples(index=False):
        nid = row.nct_id
        iname = row.name_upper
        matched_drugs = intervention_mapping.get(iname, [])
        for dname in matched_drugs:
            if nid in nct_to_cond:
                for cname in nct_to_cond[nid]:
                    drug_to_conds[dname].add(cname)
                    
    print("Generating candidate pairs...")
    benchmarks = {'SILDENAFIL', 'MINOXIDIL', 'METFORMIN', 'THALIDOMIDE', 'DULOXETINE', 'GABAPENTIN', 'BUPROPION', 'RALOXIFENE'}
    
    candidates = []
    for dname, conds in drug_to_conds.items():
        if dname in benchmarks: continue
        ind = drugs[dname]
        for cname in conds:
            if len(cname) > 4 and not is_circular(cname, ind):
                candidates.append({
                    "drug": dname,
                    "indication": ind,
                    "condition": cname
                })
                
    # Sample down to 500 for the LLM pass
    if len(candidates) > 500:
        random.seed(42)
        candidates = random.sample(candidates, 500)
        
    with open('candidates.json', 'w') as f:
        json.dump(candidates, f, indent=2)
        
    print(f"Saved {len(candidates)} candidates to candidates.json")

if __name__ == "__main__":
    main()
