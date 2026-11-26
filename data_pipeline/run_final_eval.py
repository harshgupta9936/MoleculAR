import numpy as np
import pickle
import xgboost as xgb
import pipeline
from scoring_engine import repurposing_score_bulk

VALIDATION_DRUGS = {"SILDENAFIL", "MINOXIDIL", "METFORMIN", "THALIDOMIDE",
                    "DULOXETINE", "GABAPENTIN", "BUPROPION", "RALOXIFENE"}

CASES = [
    {"drug": "SILDENAFIL", "repurposed_for": "Pulmonary Arterial Hypertension", "expected": ["TADALAFIL", "VARDENAFIL", "AVANAFIL", "RIOCIGUAT"]},
    {"drug": "MINOXIDIL", "repurposed_for": "Hair Loss", "expected": ["FINASTERIDE", "DUTASTERIDE"]},
    {"drug": "METFORMIN", "repurposed_for": "Polycystic Ovary Syndrome", "expected": ["CLOMIFENE", "LETROZOLE", "SPIRONOLACTONE"]},
    {"drug": "THALIDOMIDE", "repurposed_for": "Multiple Myeloma", "expected": ["LENALIDOMIDE", "POMALIDOMIDE", "BORTEZOMIB"]},
    {"drug": "DULOXETINE", "repurposed_for": "Diabetic Neuropathy", "expected": ["PREGABALIN", "GABAPENTIN", "AMITRIPTYLINE"]},
    {"drug": "GABAPENTIN", "repurposed_for": "Neuropathic Pain", "expected": ["PREGABALIN", "DULOXETINE", "AMITRIPTYLINE"]},
    {"drug": "BUPROPION", "repurposed_for": "Smoking Cessation", "expected": ["VARENICLINE", "CYTISINE", "NICOTINE"]},
    {"drug": "RALOXIFENE", "repurposed_for": "Osteoporosis", "expected": ["ALENDRONIC ACID", "DENOSUMAB", "TERIPARATIDE", "TAMOXIFEN"]}
]

def rank_ml_corrected(query_drug, model):
    q_feats = pipeline.drug_features.get(query_drug)
    if not q_feats:
        return []
        
    candidates = []
    for dname, feats in pipeline.drug_features.items():
        if dname == query_drug or dname in VALIDATION_DRUGS:
            continue
        # We assume pipeline.drug_features handles non-trivial ones if needed, 
        # but let's exclude trivial categories.
        if feats.get("is_trivial", False):
            continue
            
        candidates.append(feats)
        
    results = repurposing_score_bulk(q_feats, candidates, 1, 1, 1)
    
    scores = []
    X = []
    for i, c in enumerate(candidates):
        r = results[i]
        X.append([r["struct_sim"], r["target_sim"], r["side_effect_sim"]])
        
    if not X: return []
    X = np.array(X)
    probs = model.predict_proba(X)[:, 1]
    
    for i, c in enumerate(candidates):
        dname = c["drug_name"]
        pct_rank = c.get("hubness_percentile", 0.5)
        if pct_rank is None: pct_rank = 0.5
        penalty_factor = 1.0 + (pct_rank ** 2)
        corrected_prob = probs[i] / penalty_factor
        scores.append((dname, corrected_prob))
        
    scores.sort(key=lambda x: x[1], reverse=True)
    return scores[:10]

def judge_hit(ml_names, expected_list):
    for expected in expected_list:
        if expected in ml_names[:5]: return "Hit"
    for expected in expected_list:
        if expected in ml_names[:10]: return "Partial"
    return "Miss"

def main():
    print("Loading XGBoost model...")
    with open("repurposing_model.pkl", "rb") as f:
        model = pickle.load(f)
        
    print("\n" + "="*50)
    print("FINAL 8-CASE HELD-OUT VALIDATION")
    print("="*50)
    
    hits, partials, misses = 0, 0, 0
    for case in CASES:
        drug = case["drug"]
        ml_corr = rank_ml_corrected(drug, model)
        if not ml_corr:
            print(f"{drug}: Not found in DB")
            continue
            
        ml_names = [n for n, _ in ml_corr]
        judgment = judge_hit(ml_names, case["expected"])
        
        if judgment == "Hit": hits += 1
        elif judgment == "Partial": partials += 1
        else: misses += 1
        
        print(f"\nCASE: {drug} -> {case['repurposed_for']}")
        print(f"Result: {judgment}")
        print(f"Expected to find: {', '.join(case['expected'])}")
        print(f"Top 10 ML Candidates:")
        for i, (name, prob) in enumerate(ml_corr[:10], 1):
            print(f"  {i}. {name} (Prob: {prob:.4f})")
            
    print("\n" + "="*50)
    print(f"FINAL ML PERFORMANCE SUMMARY")
    print(f"Hits: {hits}/8")
    print(f"Partials: {partials}/8")
    print(f"Misses: {misses}/8")
    print("="*50)

if __name__ == "__main__":
    main()
