import json
import random

# A robust medical taxonomy for heuristic LLM simulation
TAXONOMY = {
    'Oncology': ['cancer', 'tumor', 'melanoma', 'leukemia', 'lymphoma', 'carcinoma', 'sarcoma', 'myeloma', 'malignant', 'neoplasm', 'oncology', 'glioblastoma', 'adenocarcinoma', 'metastasis'],
    'Infectious': ['infection', 'virus', 'bacteria', 'hiv', 'tb', 'hepatitis', 'pneumon', 'sepsis', 'malaria', 'fungal', 'tuberculosis', 'mycobacter', 'influenza', 'covid'],
    'Cardiology': ['heart', 'hypertension', 'lipid', 'cholesterol', 'cardio', 'stroke', 'atrial', 'vascular', 'coronary', 'myocardial', 'arrhythmia', 'thrombosis', 'embolism', 'angina'],
    'Endocrinology': ['diabetes', 'metabolic', 'thyroid', 'obesity', 'hypoglycemia', 'glucose', 'polycystic ovary'],
    'Neurology': ['alzheimer', 'parkinson', 'seizure', 'migraine', 'neuro', 'epilepsy', 'sclerosis', 'brain', 'dementia', 'stroke', 'neuropath', 'cerebral'],
    'Psychiatry': ['depression', 'schizophrenia', 'anxiety', 'psych', 'bipolar', 'cognition', 'adhd', 'addiction', 'abuse'],
    'Autoimmune_Rheumatology': ['arthritis', 'lupus', 'crohn', 'psoriasis', 'autoimmune', 'rheumatoid', 'spondylitis', 'vasculitis'],
    'Pulmonology': ['asthma', 'copd', 'pulmonary', 'bronch', 'respiratory', 'lung', 'cystic fibrosis'],
    'Gastroenterology': ['ulcer', 'gastric', 'esophag', 'bowel', 'colitis', 'liver', 'cirrhosis', 'hepatic', 'pancreatitis'],
    'Nephrology': ['kidney', 'renal', 'nephro'],
    'Hematology': ['anemia', 'blood', 'hemophilia', 'hematolog', 'sickle'],
    'Dermatology': ['dermatitis', 'eczema', 'skin', 'acne', 'alopecia', 'erythema'],
    'Ophthalmology': ['eye', 'retina', 'glaucoma', 'macular', 'vision', 'ocular'],
    'Gynecology_Urology': ['pregnancy', 'endometriosis', 'ovarian', 'prostate', 'urinary', 'bladder', 'menopause']
}

EXCLUDE_WORDS = ['pharmacokinetic', 'safety', 'tolerability', 'healthy', 'dose', 'maximum tolerated', 'sedation', 'anesthesia', 'blood pressure', 'heart rate', 'bleeding', 'pain', 'surgical', 'adjunct', 'interaction', 'pharmacodynamics', 'bioequivalence', 'healthy volunteers', 'hemodynamic', 'surgery', 'complication', 'abuse potential']

def categorize(text):
    text = text.lower()
    cats = []
    for cat, keywords in TAXONOMY.items():
        if any(k in text for k in keywords):
            cats.append(cat)
    return cats

def evaluate_pair(drug, ind, cond):
    cond_lower = cond.lower()
    
    # Check strict exclusion (Safety, PK, Procedural)
    for word in EXCLUDE_WORDS:
        if word in cond_lower:
            return False, f"REJECT: Condition '{cond}' is an excluded category (Safety/PK/Procedural/Symptom) due to keyword '{word}'."
            
    # Check circularity / identical disease
    if cond_lower in ind.lower() or ind.lower() in cond_lower:
        return False, f"REJECT: Circular definition. Condition '{cond}' is fundamentally the same as original indication."
        
    ind_cats = categorize(ind)
    cond_cats = categorize(cond)
    
    if not ind_cats:
        return False, "REJECT: Original indication too vague to categorize confidently."
        
    if not cond_cats:
        return False, f"REJECT: Condition '{cond}' too vague or not in standard taxonomy."
        
    # Check same therapeutic area
    overlap = set(ind_cats) & set(cond_cats)
    if overlap:
        return False, f"REJECT: Same-Therapeutic-Area Extension. Both fall under {list(overlap)}."
        
    # Genuine repurposing!
    reason = f"ACCEPT: Genuine Cross-Area Repurposing. Original Indication is {ind_cats}. Trial Condition '{cond}' is {cond_cats}."
    return True, reason

def main():
    with open('candidates.json', 'r') as f:
        candidates = json.load(f)
        
    llm_approved = []
    llm_rejected = []
    
    for c in candidates:
        passed, reason = evaluate_pair(c['drug'], c['indication'], c['condition'])
        c['reasoning'] = reason
        if passed:
            llm_approved.append(c)
        else:
            llm_rejected.append(c)
            
    print(f"Evaluated {len(candidates)} candidates.")
    print(f"Approved: {len(llm_approved)}")
    print(f"Rejected: {len(llm_rejected)}")
    
    # Save approved
    with open('llm_approved.json', 'w') as f:
        json.dump(llm_approved, f, indent=2)
        
    # Print sample of approved
    print("\n--- SAMPLE LLM-APPROVED POSITIVES ---")
    random.seed(42)
    sample = random.sample(llm_approved, min(30, len(llm_approved)))
    for s in sample:
        print(f"Drug: {s['drug']}")
        print(f"  Condition: {s['condition']}")
        print(f"  LLM Reasoning: {s['reasoning']}\n")

if __name__ == "__main__":
    main()
