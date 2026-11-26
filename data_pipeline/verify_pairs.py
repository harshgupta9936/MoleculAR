import json
import sqlite3
import re

def get_broad_category(text):
    text = text.lower()
    if any(w in text for w in ['cancer', 'tumor', 'melanoma', 'leukemia', 'lymphoma', 'carcinoma', 'sarcoma', 'myeloma', 'malignant', 'neoplasm']):
        return 'Oncology'
    if any(w in text for w in ['infection', 'virus', 'bacteria', 'hiv', 'tb', 'hepatitis', 'fung', 'mycobacter', 'pneumonia', 'sepsis']):
        return 'Infectious Disease'
    if any(w in text for w in ['heart', 'hypertension', 'diabetes', 'lipid', 'cholesterol', 'cardio', 'stroke', 'atrial', 'vascular']):
        return 'Cardiovascular/Metabolic'
    if any(w in text for w in ['depression', 'schizophrenia', 'alzheimer', 'parkinson', 'seizure', 'pain', 'migraine', 'anxiety', 'psych', 'neuro']):
        return 'Neurological/Psychiatric'
    if any(w in text for w in ['arthritis', 'lupus', 'crohn', 'psoriasis', 'asthma', 'allergy', 'autoimmune', 'inflam']):
        return 'Autoimmune/Inflammatory'
    return 'Other/Unknown'

def main():
    conn = sqlite3.connect('drugrepo.db')
    c = conn.cursor()
    
    with open('training_pairs.json', 'r') as f:
        data = json.load(f)
        
    positives = data['positives']
    negatives = data['negatives']
    
    # 1. Classify Positives
    same_area = 0
    cross_area = 0
    
    for dname, cond in positives:
        c.execute("SELECT indication FROM drugs WHERE generic_base_name = ? LIMIT 1", (dname,))
        res = c.fetchone()
        ind = res[0] if res and res[0] else ''
        
        cat1 = get_broad_category(ind)
        cat2 = get_broad_category(cond)
        
        # Word overlap check for very direct extensions (e.g., "type 1 diabetes" vs "type 2 diabetes")
        words1 = set(re.findall(r'\b\w{5,}\b', ind.lower()))
        words2 = set(re.findall(r'\b\w{5,}\b', cond.lower()))
        
        if (cat1 != 'Other/Unknown' and cat1 == cat2) or (words1 & words2):
            same_area += 1
        else:
            cross_area += 1
            
    print(f"Total Positives: {len(positives)}")
    print(f"Same-Therapeutic-Area Extension: {same_area} ({same_area/len(positives)*100:.1f}%)")
    print(f"Cross-Therapeutic-Area: {cross_area} ({cross_area/len(positives)*100:.1f}%)")
    
    # 2. Deep Dive on DIRITHROMYCIN
    print("\n--- DIRITHROMYCIN DEEP DIVE ---")
    dname = "DIRITHROMYCIN"
    cond = "tuberculosis, pulmonary, drug sensitive"
    
    # Get targets of DIRITHROMYCIN
    c.execute("SELECT target_name FROM targets t JOIN drugs d ON t.chembl_id = d.chembl_id WHERE d.generic_base_name = ?", (dname,))
    diri_targets = set([r[0].lower() for r in c.fetchall()])
    print(f"DIRITHROMYCIN targets: {diri_targets}")
    
    # In my build script, I got tested drugs from AACT mapping. Let's recreate tested drugs for this condition.
    # To do this correctly, we simulate what build_training_data did:
    import pandas as pd
    cond_df = pd.read_csv('conditions.txt', sep='|', usecols=['nct_id', 'downcase_name'], on_bad_lines='skip', low_memory=False)
    nct_ids_for_tb = cond_df[cond_df['downcase_name'] == cond]['nct_id'].tolist()
    print(f"Trials for '{cond}': {len(nct_ids_for_tb)}")
    
    int_df = pd.read_csv('interventions.txt', sep='|', usecols=['nct_id', 'intervention_type', 'name'], on_bad_lines='skip', low_memory=False)
    tb_interventions = int_df[(int_df['nct_id'].isin(nct_ids_for_tb)) & (int_df['intervention_type'].str.upper() == 'DRUG')]['name'].tolist()
    
    print(f"Raw interventions tested for this TB condition: {list(set(tb_interventions))[:10]}")
    
    # Which of our active drugs map to these interventions?
    c.execute("SELECT generic_base_name FROM drugs WHERE trivial_category IS NULL")
    active_drugs = [r[0] for r in c.fetchall()]
    
    tb_tested_drugs = set()
    for iname in tb_interventions:
        iname_upper = str(iname).upper()
        for active in active_drugs:
            if len(active) >= 5 and active in iname_upper:
                tb_tested_drugs.add(active)
                
    print(f"Active dataset drugs mapped to this condition: {tb_tested_drugs}")
    
    # Get all targets for these tb_tested_drugs
    tb_targets = set()
    for tbd in tb_tested_drugs:
        c.execute("SELECT target_name FROM targets t JOIN drugs d ON t.chembl_id = d.chembl_id WHERE d.generic_base_name = ?", (tbd,))
        for r in c.fetchall(): tb_targets.add(r[0].lower())
        
    print(f"Total distinct targets of these TB drugs: {len(tb_targets)}")
    
    overlap = diri_targets & tb_targets
    print(f"Target overlap: {overlap}")
    
    conn.close()

if __name__ == "__main__":
    main()
