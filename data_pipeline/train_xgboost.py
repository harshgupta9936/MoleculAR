import json
import csv
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
import xgboost as xgb
import pipeline
from scoring_engine import repurposing_score_bulk
import pickle

def load_data():
    with open("disease_references.json", "r") as f:
        valid_diseases = json.load(f)
        
    df = pd.read_csv("training_labels.csv")
    
    X = []
    y = []
    meta = []
    
    print("Computing feature vectors for dataset...")
    for idx, row in df.iterrows():
        drug = row["drug_name"]
        disease = row["disease"]
        label = row["label"]
        
        q_feats = pipeline.drug_features.get(drug)
        if not q_feats:
            continue
            
        refs = valid_diseases.get(disease, [])
        r_feats_list = [pipeline.drug_features.get(r) for r in refs if pipeline.drug_features.get(r)]
        if not r_feats_list:
            continue
            
        bulk_results = repurposing_score_bulk(q_feats, r_feats_list, w_structure=0.4, w_target=0.4, w_sideeffect=0.2)
        
        max_struct = max([r["struct_sim"] for r in bulk_results])
        max_target = max([r["target_sim"] for r in bulk_results])
        max_se = max([r["side_effect_sim"] for r in bulk_results])
        
        X.append([max_struct, max_target, max_se])
        y.append(label)
        meta.append((drug, disease))
        
    return np.array(X), np.array(y), meta

def train_and_evaluate():
    X, y, meta = load_data()
    
    print(f"\nDataset loaded: {len(X)} total pairs ({sum(y)} positive, {len(y)-sum(y)} negative)")
    print("\nStarting Stratified 5-Fold Cross-Validation...")
    
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    fold_aucs = []
    oof_predictions = {}
    importances = []
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y), 1):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        model = xgb.XGBClassifier(
            n_estimators=100, 
            max_depth=3, 
            learning_rate=0.1, 
            random_state=42,
            eval_metric='logloss'
        )
        model.fit(X_train, y_train)
        
        preds = model.predict_proba(X_test)[:, 1]
        auc = roc_auc_score(y_test, preds)
        fold_aucs.append(auc)
        importances.append(model.feature_importances_)
        
        print(f"  Fold {fold} AUC: {auc:.4f}")
        
        for i, idx in enumerate(test_idx):
            oof_predictions[meta[idx]] = preds[i]
            
    print(f"\nMean ROC-AUC: {np.mean(fold_aucs):.4f} (Std: {np.std(fold_aucs):.4f})")
    
    mean_imp = np.mean(importances, axis=0)
    print("\nFeature Importances:")
    print(f"  struct_sim:      {mean_imp[0]:.4f}")
    print(f"  target_sim:      {mean_imp[1]:.4f}")
    print(f"  side_effect_sim: {mean_imp[2]:.4f}")
    
    print("\nEvaluation Cases (Out-of-Fold Predicted Probabilities):")
    eval_cases = [
        ("SILDENAFIL", "Pulmonary Arterial Hypertension"),
        ("GABAPENTIN", "Neuropathic Pain"),
        ("DULOXETINE", "Neuropathic Pain")
    ]
    
    for drug, disease in eval_cases:
        if (drug, disease) in oof_predictions:
            prob = oof_predictions[(drug, disease)]
            print(f"  {drug} -> {disease}: {prob:.4f}")
        else:
            print(f"  {drug} -> {disease}: Not found in dataset")
            
    print("\nTraining final model on all data and saving...")
    final_model = xgb.XGBClassifier(
        n_estimators=100, 
        max_depth=3, 
        learning_rate=0.1, 
        random_state=42,
        eval_metric='logloss'
    )
    final_model.fit(X, y)
    with open("repurposing_model.pkl", "wb") as f:
        pickle.dump(final_model, f)
    print("Saved repurposing_model.pkl")

if __name__ == "__main__":
    # Ensure SQLite DB is initialized inside pipeline before running
    train_and_evaluate()
