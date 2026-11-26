import pandas as pd
import numpy as np
import json
import pickle
import xgboost as xgb
from scoring_engine import get_fingerprint, repurposing_score
from score_hit_rate import judge_hit, compute_diversity, CASES

VALIDATION_DRUGS = {"SILDENAFIL", "MINOXIDIL", "METFORMIN", "THALIDOMIDE",
                    "DULOXETINE", "GABAPENTIN", "BUPROPION", "RALOXIFENE"}

def load_data():
    df = pd.read_csv("full_dataset.csv")
    with open("drug_targets.json", "r") as f:
        targets = json.load(f)
    with open("drug_side_effects.json", "r") as f:
        side_effects = json.load(f)
        
    def extract_targets(t_list):
        if not t_list: return set()
        if isinstance(t_list[0], dict): return set(t["target_id"] for t in t_list)
        return set(t_list)

    drug_features = {}
    for _, row in df.iterrows():
        drug = row["drug_name"]
        drug_features[drug] = {
            "fingerprint": get_fingerprint(row["smiles"]),
            "targets": extract_targets(targets.get(row["chembl_id"], [])),
            "side_effects": set(side_effects.get(drug, []))
        }
    return df, drug_features

def compute_hubness(df, drug_features):
    avg_sim = {}
    for drug in df["drug_name"]:
        if drug not in drug_features: continue
        feats1 = drug_features[drug]
        scores = []
        for other in df["drug_name"]:
            if other != drug and other in drug_features:
                feats2 = drug_features[other]
                s = repurposing_score(feats1, feats2)["total_score"]
                scores.append(s)
        avg_sim[drug] = np.mean(scores)
    
    hub_pct_ranks = pd.Series(avg_sim).rank(pct=True)
    return hub_pct_ranks

def rank_ml_corrected(query_drug, drug_features, df, hub_pct_ranks, model):
    feats1 = drug_features[query_drug]
    scores = []
    for other in df["drug_name"]:
        if other == query_drug or other not in drug_features or other in VALIDATION_DRUGS:
            continue
        feats2 = drug_features[other]
        
        sims = repurposing_score(feats1, feats2, 1, 1, 1)
        X = np.array([[sims["struct_sim"], sims["target_sim"], sims["side_effect_sim"]]])
        
        raw_prob = model.predict_proba(X)[0][1]
        
        pct_rank = hub_pct_ranks.get(other, 0.5)
        penalty_factor = 1.0 + (pct_rank ** 2)
        corrected_prob = raw_prob / penalty_factor
        
        scores.append((other, corrected_prob))
        
    scores.sort(key=lambda x: x[1], reverse=True)
    return scores[:10]

def main():
    print("Loading data...")
    df, drug_features = load_data()
    print("Computing hubness percentiles for correction...")
    hub_pct_ranks = compute_hubness(df, drug_features)
    
    with open("repurposing_model.pkl", "rb") as f:
        model = pickle.load(f)
        
    print("Scoring 8 cases using ML model...\n")
    rows = []
    
    for case in CASES:
        drug = case["drug"]
        if drug not in drug_features: continue
        
        ml_corr = rank_ml_corrected(drug, drug_features, df, hub_pct_ranks, model)
        ml_names = [n for n, _ in ml_corr]
        
        judgment, reason = judge_hit(drug, ml_names, case)
        diversity = compute_diversity(ml_corr, drug_features)
        
        rows.append({
            "case": f"{drug} -> {case['repurposed_for']}",
            "result": judgment,
            "reason": reason,
            "diversity": diversity,
            "top3": ", ".join(ml_names[:3])
        })
        
    # Summarize Hand-Weighted Phase 4.4 Results for comparison
    hw_hits, hw_partials, hw_misses = 3, 4, 1
    hw_div = 0.1269
    
    ml_hits = sum(1 for r in rows if r['result'] == 'hit')
    ml_partials = sum(1 for r in rows if r['result'] == 'partial')
    ml_misses = sum(1 for r in rows if r['result'] == 'miss')
    ml_div = np.mean([r['diversity'] for r in rows])
    
    with open("ml_vs_handweighted_comparison.md", "w", encoding="utf-8") as out:
        out.write("# ML Classifier vs Hand-Weighted Model\n\n")
        out.write("This compares the hit rate of the XGBoost ML model (trained on the filtered 22-row dataset) ")
        out.write("against the Phase 4.4 Hand-Weighted model (0.4/0.4/0.2 weights). Both models use the ")
        out.write("percentile-rank hub bias correction.\n\n")
        
        out.write("## Summary Comparison\n\n")
        out.write("| Metric | Hand-Weighted (Phase 4.4) | ML XGBoost (Phase 6) |\n")
        out.write("|---|---|---|\n")
        out.write(f"| Hits | {hw_hits}/8 ({hw_hits/8*100:.0f}%) | {ml_hits}/8 ({ml_hits/8*100:.0f}%) |\n")
        out.write(f"| Partials | {hw_partials}/8 ({hw_partials/8*100:.0f}%) | {ml_partials}/8 ({ml_partials/8*100:.0f}%) |\n")
        out.write(f"| Misses | {hw_misses}/8 ({hw_misses/8*100:.0f}%) | {ml_misses}/8 ({ml_misses/8*100:.0f}%) |\n")
        out.write(f"| Hit+Partial rate | {(hw_hits+hw_partials)/8*100:.0f}% | {(ml_hits+ml_partials)/8*100:.0f}% |\n")
        out.write(f"| Avg intra-list diversity | {hw_div:.4f} | {ml_div:.4f} |\n\n")
        
        out.write("## ML Per-Case Detail\n\n")
        for r in rows:
            out.write(f"### {r['case']}\n")
            out.write(f"- **Result ({r['result']}):** {r['reason']}\n")
            out.write(f"- **Top 3 Candidates:** {r['top3']}\n\n")
            
    print("Saved ml_vs_handweighted_comparison.md")

if __name__ == "__main__":
    main()
