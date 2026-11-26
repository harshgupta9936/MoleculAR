import pandas as pd
import sqlite3
import random
import re
import json
import itertools

TAXONOMY = {
    'Oncology': ['cancer', 'tumor', 'melanoma', 'leukemia', 'lymphoma', 'carcinoma', 'sarcoma', 'myeloma', 'malignant', 'neoplasm', 'oncology', 'glioblastoma', 'adenocarcinoma', 'metastasis'],
    'Infectious': ['infection', 'virus', 'bacteria', 'hiv', 'tb', 'hepatitis', 'pneumon', 'sepsis', 'malaria', 'fungal', 'tuberculosis', 'mycobacter', 'influenza', 'covid', 'bacteremia', 'rsv', 'viral'],
    'Cardiology': ['heart', 'hypertension', 'lipid', 'cholesterol', 'cardio', 'stroke', 'atrial', 'vascular', 'coronary', 'myocardial', 'arrhythmia', 'thrombosis', 'embolism', 'angina', 'fibrillation', 'failure', 'infarction'],
    'Endocrinology': ['diabetes', 'metabolic', 'thyroid', 'obesity', 'hypoglycemia', 'hyperglycemia', 'blood sugar', 'glucose', 'polycystic ovary', 'endocrine'],
    'Neurology': ['alzheimer', 'parkinson', 'seizure', 'migraine', 'neuro', 'epilepsy', 'sclerosis', 'brain', 'dementia', 'stroke', 'neuropath', 'cerebral', 'cognitive'],
    'Psychiatry': ['depression', 'schizophrenia', 'anxiety', 'psych', 'bipolar', 'cognition', 'adhd', 'addiction', 'abuse', 'disorder', 'panic', 'insomnia'],
    'Autoimmune_Rheumatology': ['arthritis', 'lupus', 'crohn', 'psoriasis', 'autoimmune', 'rheumatoid', 'spondylitis', 'vasculitis', 'inflammatory', 'erythematosus'],
    'Pulmonology': ['asthma', 'copd', 'pulmonary', 'bronch', 'respiratory', 'lung', 'cystic fibrosis', 'fibrosis'],
    'Gastroenterology': ['ulcer', 'gastric', 'esophag', 'bowel', 'colitis', 'liver', 'cirrhosis', 'hepatic', 'pancreatitis', 'gastrointestinal'],
    'Nephrology': ['kidney', 'renal', 'nephro', 'nephritis', 'glomerulonephritis'],
    'Hematology': ['anemia', 'blood', 'hemophilia', 'hematolog', 'sickle', 'thrombocytopenia', 'coagulation', 'engraftment'],
    'Dermatology': ['dermatitis', 'eczema', 'skin', 'acne', 'alopecia', 'erythema', 'psoriasis'],
    'Ophthalmology': ['eye', 'retina', 'glaucoma', 'macular', 'vision', 'ocular', 'uveitis'],
    'Gynecology_Urology': ['pregnancy', 'endometriosis', 'ovarian', 'prostate', 'urinary', 'bladder', 'menopause', 'cervical']
}

def categorize(text):
    text = text.lower()
    cats = []
    for cat, keywords in TAXONOMY.items():
        if any(re.search(rf'\b{k}\b', text) for k in keywords):
            cats.append(cat)
    return cats

def is_valid_drug_match(iname, dname):
    return bool(re.search(rf'\b{re.escape(dname)}\b', str(iname), re.IGNORECASE))

def main():
    print("Loading valid positives...")
    with open('heuristic_approved_v8.json') as f:
        positives = json.load(f)
        
    num_pos = len(positives)
    target_neg = num_pos * 2
    print(f"Targeting {target_neg} negatives to match 1:2 ratio with {num_pos} positives.")
    
    print("Loading drugs...")
    conn = sqlite3.connect('drugrepo.db')
    c = conn.cursor()
    c.execute("SELECT chembl_id, generic_base_name, indication FROM drugs WHERE indication != 'Other' AND indication IS NOT NULL AND indication != '' AND trivial_category IS NULL")
    
    drugs = {}
    for cid, name, ind in c.fetchall():
        dname = str(name).strip().upper()
        if len(dname) >= 4:
            drugs[dname] = str(ind).lower()
    conn.close()
    drug_names = list(drugs.keys())
    
    print("Loading conditions and interventions...")
    # Load ALL conditions from AACT, not just treatment phase ones, because we want to completely exclude ANY off-label/investigational evidence.
    cond_df = pd.read_csv('conditions.txt', sep='|', usecols=['nct_id', 'downcase_name'], on_bad_lines='skip', low_memory=False)
    int_df = pd.read_csv('interventions.txt', sep='|', usecols=['nct_id', 'intervention_type', 'name'], on_bad_lines='skip', low_memory=False)
    int_df = int_df[int_df['intervention_type'].astype(str).str.upper().isin(['DRUG', 'BIOLOGICAL', 'DIETARY SUPPLEMENT'])]
    
    print("Mapping drugs to NCTs (global AACT coverage)...")
    unique_interventions = int_df['name'].dropna().unique()
    intervention_mapping = {}
    for iname in unique_interventions:
        matched = []
        for dname in drug_names:
            if is_valid_drug_match(iname, dname):
                matched.append(dname)
        if matched:
            intervention_mapping[iname] = matched
            
    drug_to_ncts = {d: set() for d in drug_names}
    for row in int_df.itertuples(index=False):
        nid = row.nct_id
        iname = row.name
        for dname in intervention_mapping.get(iname, []):
            drug_to_ncts[dname].add(nid)
            
    nct_to_cond = {}
    cond_counts = {}
    for row in cond_df.itertuples(index=False):
        nid = row.nct_id
        cname = str(row.downcase_name).lower()
        if cname == 'nan': continue
        if nid not in nct_to_cond: nct_to_cond[nid] = set()
        nct_to_cond[nid].add(cname)
        cond_counts[cname] = cond_counts.get(cname, 0) + 1
        
    print("Building set of existing drug-condition pairs from ANY trial...")
    aact_pairs = set()
    for dname, ncts in drug_to_ncts.items():
        for nid in ncts:
            if nid in nct_to_cond:
                for cname in nct_to_cond[nid]:
                    if len(cname) > 4:
                        aact_pairs.add((dname, cname))
                        
    legacy_df = pd.read_csv('training_labels.csv', dtype=str, on_bad_lines='skip')
    legacy_pairs = set()
    for _, row in legacy_df.iterrows():
        d = str(row['drug_name']).upper()
        c = str(row['disease']).lower()
        legacy_pairs.add((d, c))
        
    benchmarks = {'SILDENAFIL', 'MINOXIDIL', 'METFORMIN', 'THALIDOMIDE', 'DULOXETINE', 'GABAPENTIN', 'BUPROPION', 'RALOXIFENE'}
    
    # We want a pool of candidate conditions that are robust (not trivially all-zero feature vectors)
    # i.e., conditions that appear in at least 10 distinct trials (so they have some drug associations)
    robust_conditions = [c for c, count in cond_counts.items() if count >= 10 and len(c) > 4]
    print(f"Robust condition pool size (>=10 trials): {len(robust_conditions)}")
    
    negatives = []
    attempts = 0
    
    HOT_CLASSES = ['STATIN', 'FIBRATE', 'GLIFLOZIN', 'METFORMIN', 'ASPIRIN', 'CELECOXIB', 'SILDENAFIL', 'TADALAFIL']
    HOT_TARGETS = ['Oncology', 'Neurology', 'Infectious']
    
    print("Generating negatives...")
    random.seed(99) # For reproducibility
    
    while len(negatives) < target_neg and attempts < 1000000:
        attempts += 1
        dname = random.choice(drug_names)
        if dname in benchmarks: continue
        
        cname = random.choice(robust_conditions)
        
        # 1. Reject if pair was in AACT
        if (dname, cname) in aact_pairs: continue
        
        # 2. Reject if in legacy pairs
        if (dname, cname) in legacy_pairs: continue
        
        ind = drugs[dname]
        
        # 3. String overlap checks
        if cname in ind or ind in cname: continue
        words1 = set(re.findall(r'\b\w{5,}\b', cname))
        words2 = set(re.findall(r'\b\w{5,}\b', ind))
        if words1 & words2: continue
        
        # 4. Taxonomy check
        ind_cats = categorize(ind)
        cond_cats = categorize(cname)
        
        if not ind_cats or not cond_cats: continue
        if set(ind_cats) & set(cond_cats): continue
            
        # 5. Prevent literature repurposing leak (Statin/Fibrate/Gliflozin to Cancer/Neuro/Infectious)
        is_hot_drug = any(h in dname for h in HOT_CLASSES)
        is_hot_target = any(t in cond_cats for t in HOT_TARGETS)
        if is_hot_drug and is_hot_target: continue
        
        negatives.append({
            "drug": dname,
            "condition": cname,
            "indication": ind,
            "reasoning": f"Original: {ind_cats} | Cond: {cond_cats}"
        })
        
    print(f"Generated {len(negatives)} negative pairs in {attempts} attempts.")
    
    with open('generated_negatives.json', 'w') as f:
        json.dump(negatives, f, indent=2)
        
    print("\n==================================================")
    print("20-SAMPLE NEGATIVES AUDIT (RAW OUTPUT)")
    print("==================================================")
    
    for i, s in enumerate(negatives[:20]):
        print(f"\n{i+1}. Pair: {s['drug']} -> {s['condition']}")
        print(f"   Indication: {s['indication']}")
        print(f"   Categories: {s['reasoning']}")
        print(f"   Feature Sanity: Condition '{s['condition']}' appears in {cond_counts[s['condition']]} AACT trials.")

if __name__ == "__main__":
    main()
