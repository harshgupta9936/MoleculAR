import pandas as pd
import random

# Seed for reproducibility in negative sampling
random.seed(42)

# Strict FDA-Approved Multi-Indication Drugs ONLY
# Format: Drug -> [List of genuinely distinct FDA-approved indications]
# Source: FDA labels, manually verified.
# NOTE: Metformin is excluded per rules (only one accepted FDA indication in this context).
POSITIVE_CASES = {
    # The 7 valid multi-indication validation cases
    "SILDENAFIL": ["Erectile Dysfunction", "Pulmonary Arterial Hypertension"],
    "MINOXIDIL": ["Hypertension", "Hair Loss"],
    "THALIDOMIDE": ["Multiple Myeloma", "Leprosy (ENL)"],
    "DULOXETINE": ["Depression", "Diabetic Neuropathy"],
    "GABAPENTIN": ["Epilepsy", "Neuropathic Pain (PHN)"],
    "BUPROPION": ["Depression", "Smoking Cessation"],
    "RALOXIFENE": ["Osteoporosis", "Breast Cancer Prevention"],
    
    # 4 Additional manually curated multi-indication drugs present in the 183 dataset
    "VALPROIC ACID": ["Epilepsy", "Bipolar Mania"],
    "CARBAMAZEPINE": ["Epilepsy", "Trigeminal Neuralgia"],
    "TERAZOSIN": ["Hypertension", "Benign Prostatic Hyperplasia"],
    "FINASTERIDE": ["Benign Prostatic Hyperplasia", "Hair Loss"]
}

def build_dataset():
    # Load dataset to verify drugs exist
    df = pd.read_csv("full_dataset.csv")
    valid_drugs = set(df["drug_name"].str.upper())
    
    records = []
    all_indications = set()
    
    # 1. Build Positive Class
    for drug, indications in POSITIVE_CASES.items():
        if drug not in valid_drugs:
            print(f"WARNING: {drug} not in dataset. Skipping.")
            continue
            
        for ind in indications:
            records.append({
                "drug_name": drug,
                "disease": ind,
                "label": 1,
                "source": "FDA Approved Indication"
            })
            all_indications.add(ind)
            
    # 2. Build Negative Class
    # We want a balanced dataset (1:1 positive:negative ratio)
    # We will assign random indications to these same multi-indication drugs,
    # ensuring they are NOT actually approved for that indication.
    positive_count = len(records)
    all_indications_list = list(all_indications)
    
    drugs_used = list(POSITIVE_CASES.keys())
    negative_records = []
    
    while len(negative_records) < positive_count:
        drug = random.choice(drugs_used)
        rand_ind = random.choice(all_indications_list)
        
        # Ensure it's a true negative (not in its known approved list)
        if rand_ind not in POSITIVE_CASES[drug]:
            # Also ensure we haven't already added this negative pair
            pair_exists = any(r["drug_name"] == drug and r["disease"] == rand_ind for r in negative_records)
            if not pair_exists:
                negative_records.append({
                    "drug_name": drug,
                    "disease": rand_ind,
                    "label": 0,
                    "source": "Random Negative Pairing"
                })
                
    records.extend(negative_records)
    
    # Shuffle the dataset
    random.shuffle(records)
    
    # Convert to DataFrame and save
    out_df = pd.DataFrame(records)
    out_df.to_csv("training_labels.csv", index=False)
    
    print("=== Phase 5 Dataset Generation Complete ===")
    print(f"Positive pairs: {positive_count}")
    print(f"Negative pairs: {len(negative_records)}")
    print(f"Total pairs: {len(records)}")
    print("\nCAVEAT: Negative pairs are randomly matched. There is a small risk of false negatives ")
    print("(unrecognized off-label use). This is a known limitation of the methodology.")
    print("Saved to training_labels.csv")

if __name__ == "__main__":
    build_dataset()
