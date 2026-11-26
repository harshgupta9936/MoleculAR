import json
import sqlite3
import pandas as pd

TAXONOMY = {
    'Oncology': ['cancer', 'tumor', 'melanoma', 'leukemia', 'lymphoma', 'carcinoma', 'sarcoma', 'myeloma', 'malignant', 'neoplasm', 'oncology', 'glioblastoma', 'adenocarcinoma', 'metastasis', 'adenoma'],
    'Infectious': ['infection', 'virus', 'bacteria', 'hiv', 'tb', 'hepatitis', 'pneumon', 'sepsis', 'malaria', 'fungal', 'tuberculosis', 'mycobacter', 'influenza', 'covid', 'bacteremia', 'rsv', 'viral', 'antibiotic', 'antiviral'],
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
        if any(k in text for k in keywords):
            cats.append(cat)
    return cats

def main():
    conn = sqlite3.connect('drugrepo.db')
    c = conn.cursor()
    c.execute("SELECT chembl_id, generic_base_name, indication FROM drugs WHERE indication != 'Other' AND indication IS NOT NULL AND trivial_category IS NULL")
    drugs_raw = c.fetchall()
    
    mismatches = []
    
    with open('heuristic_approved.json', 'r') as f:
        approved = json.load(f)
        
    # We want to find drugs where their AACT condition is highly synonymous with their real indication, 
    # but they were mapped to different taxonomy categories.
    
    # Actually, we can just look at drugs with >= 5 trials in the approved list,
    # and print their indications and conditions to spot obvious circularity.
    
    drug_counts = {}
    drug_to_conds = {}
    drug_to_ind = {}
    
    for row in drugs_raw:
        drug_to_ind[str(row[1]).upper()] = str(row[2])
        
    for p in approved:
        d = p['drug']
        c = p['condition']
        drug_counts[d] = drug_counts.get(d, 0) + 1
        if d not in drug_to_conds: drug_to_conds[d] = set()
        drug_to_conds[d].add(c)
        
    print("Top 20 drugs by 'repurposing' trial count:")
    sorted_drugs = sorted(drug_counts.items(), key=lambda x: x[1], reverse=True)
    
    for d, count in sorted_drugs[:20]:
        print(f"\n{d} ({count} trials)")
        print(f"  DB Indication: {drug_to_ind.get(d)}")
        print(f"  Top Conditions: {list(drug_to_conds[d])[:5]}")

if __name__ == "__main__":
    main()
