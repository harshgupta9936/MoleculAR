import pandas as pd
import sqlite3
import random
import re

def load_data():
    conn = sqlite3.connect('drugrepo.db')
    cursor = conn.cursor()
    cursor.execute("SELECT chembl_id, generic_base_name, indication FROM drugs WHERE trivial_category IS NULL")
    drugs_raw = cursor.fetchall()
    
    drugs = {}
    for cid, name, ind in drugs_raw:
        dname = str(name).strip().upper()
        drugs[dname] = {
            'chembl_id': cid,
            'indication': str(ind).lower() if ind else '',
            'targets': set(),
            'side_effects': set()
        }
        
    cursor.execute("SELECT d.generic_base_name, t.target_name FROM targets t JOIN drugs d ON t.chembl_id = d.chembl_id")
    for name, tname in cursor.fetchall():
        dname = str(name).strip().upper()
        if dname in drugs: drugs[dname]['targets'].add(tname.lower())
        
    cursor.execute("SELECT d.generic_base_name, s.side_effect_name FROM side_effects s JOIN drugs d ON s.chembl_id = d.chembl_id")
    for name, se in cursor.fetchall():
        dname = str(name).strip().upper()
        if dname in drugs: drugs[dname]['side_effects'].add(se.lower())
        
    conn.close()
    return drugs

def parse_aact(drugs):
    cond_df = pd.read_csv('conditions.txt', sep='|', usecols=['nct_id', 'downcase_name'], on_bad_lines='skip', low_memory=False)
    nct_to_cond = {}
    for row in cond_df.itertuples(index=False):
        nid = row.nct_id
        cname = str(row.downcase_name).lower()
        if cname == 'nan': continue
        if nid not in nct_to_cond: nct_to_cond[nid] = set()
        nct_to_cond[nid].add(cname)
        
    int_df = pd.read_csv('interventions.txt', sep='|', usecols=['nct_id', 'intervention_type', 'name'], on_bad_lines='skip', low_memory=False)
    int_df = int_df[int_df['intervention_type'].astype(str).str.upper().isin(['DRUG', 'BIOLOGICAL', 'DIETARY SUPPLEMENT'])]
    
    drug_to_nct = {dname: set() for dname in drugs.keys()}
    nct_to_drugs = {}
    
    drug_names = list(drugs.keys())
    int_df['name_upper'] = int_df['name'].astype(str).str.upper()
    unique_interventions = int_df['name_upper'].unique()
    
    intervention_mapping = {}
    drug_names.sort(key=len, reverse=True)
    
    match_count = 0
    for iname in unique_interventions:
        if not isinstance(iname, str): continue
        matched = []
        for dname in drug_names:
            if len(dname) < 5: continue
            if dname in iname:
                matched.append(dname)
        if matched:
            match_count += 1
            intervention_mapping[iname] = matched
            
    print(f"Matched {match_count} unique interventions to our drugs.")
    
    for row in int_df.itertuples(index=False):
        nid = row.nct_id
        iname = row.name_upper
        matched_drugs = intervention_mapping.get(iname, [])
        for dname in matched_drugs:
            drug_to_nct[dname].add(nid)
            if nid not in nct_to_drugs: nct_to_drugs[nid] = set()
            nct_to_drugs[nid].add(dname)
                
    return nct_to_cond, drug_to_nct, nct_to_drugs

def is_circular(cond, ind):
    cond = cond.lower()
    ind = ind.lower()
    if cond in ind: return True
    
    # Common medical synonyms
    synonyms = [
        ('cancer', 'neoplasm'), ('cancer', 'carcinoma'), ('cancer', 'tumor'), ('cancer', 'malignan'), ('cancer', 'sarcoma'),
        ('heart', 'cardiac'), ('heart', 'myocardial'),
        ('kidney', 'renal'), ('liver', 'hepatic'), ('lung', 'pulmonary'),
        ('brain', 'cerebral'), ('stomach', 'gastric'), ('blood', 'hematologic')
    ]
    
    # Normalize condition with synonyms for checking
    for s1, s2 in synonyms:
        if s1 in cond and s2 in ind: return True
        if s2 in cond and s1 in ind: return True
        
    words1 = set(re.findall(r'\b\w{5,}\b', cond))
    words2 = set(re.findall(r'\b\w{5,}\b', ind))
    if words1 & words2: return True
    
    return False

def main():
    drugs = load_data()
    print(f"Loaded {len(drugs)} active scored drugs.")
    nct_to_cond, drug_to_nct, nct_to_drugs = parse_aact(drugs)
    
    cond_to_tested_drugs = {}
    for nid, dnames in nct_to_drugs.items():
        if nid in nct_to_cond:
            for cname in nct_to_cond[nid]:
                if cname not in cond_to_tested_drugs: cond_to_tested_drugs[cname] = set()
                for d in dnames:
                    cond_to_tested_drugs[cname].add(d)
                    
    print("Generating positive pairs...")
    positives = set()
    for dname, ncts in drug_to_nct.items():
        d_ind = drugs[dname]['indication']
        if d_ind == 'other' or not d_ind: continue # MUST have real indication to generate positives!
        
        for nid in ncts:
            if nid in nct_to_cond:
                for cname in nct_to_cond[nid]:
                    if len(cname) > 4 and not is_circular(cname, d_ind):
                        positives.add((dname, cname))
                        
    benchmarks = {'SILDENAFIL', 'MINOXIDIL', 'METFORMIN', 'THALIDOMIDE', 'DULOXETINE', 'GABAPENTIN', 'BUPROPION', 'RALOXIFENE'}
    positives = [(d, c) for (d, c) in positives if d not in benchmarks]
    if not positives:
        print("ERROR: 0 positive pairs found!")
        return
        
    positives = random.sample(list(positives), min(300, len(positives)))
    
    print("Generating negative pairs...")
    negatives = set()
    all_conditions = list(cond_to_tested_drugs.keys())
    drug_list = [d for d in drugs.keys() if d not in benchmarks and drugs[d]['indication'] != 'other' and drugs[d]['indication']]
    
    attempts = 0
    while len(negatives) < len(positives) * 1.5 and attempts < 200000:
        attempts += 1
        dname = random.choice(drug_list)
        cname = random.choice(all_conditions)
        if len(cname) <= 4: continue
        if (dname, cname) in positives: continue
        
        tested_drugs_for_c = cond_to_tested_drugs[cname]
        if dname in tested_drugs_for_c: continue
        
        if is_circular(cname, drugs[dname]['indication']): continue
        in_se = any(cname in se for se in drugs[dname]['side_effects'])
        if in_se: continue
        
        d_targets = drugs[dname]['targets']
        shared_target_found = False
        if d_targets:
            for other_d in tested_drugs_for_c:
                other_targets = drugs[other_d]['targets']
                if d_targets & other_targets:
                    shared_target_found = True
                    break
        if shared_target_found: continue
        
        negatives.add((dname, cname))
        
    negatives = list(negatives)
    print(f"\nGenerated {len(positives)} positive pairs and {len(negatives)} negative pairs.")
    
    import json
    with open('training_pairs.json', 'w') as f:
        json.dump({'positives': positives, 'negatives': negatives}, f, indent=2)
        
    overlap = [d for p in positives+negatives for d in p if d in benchmarks]
    print(f"Benchmark overlap: {len(overlap)}")
    
    print("\nExample Positives:")
    for p in positives[:5]: print(f"  {p[0]} -> {p[1]}")
    
    print("\nExample Negatives:")
    for n in negatives[:5]: print(f"  {n[0]} -> {n[1]}")

if __name__ == "__main__":
    main()
