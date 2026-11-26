import json
import pandas as pd
import sqlite3

def build():
    # 1. 9 cases from post-indication-fix audit
    post_fix = [
        {"drug": "CARMUSTINE", "condition": "autologous transplant autoimmune", "nct_id": "NCT00716066"},
        {"drug": "CHLOROQUINE", "condition": "pneumonia", "nct_id": "NCT04779047"},
        {"drug": "LAMOTRIGINE", "condition": "pathologic skin picking", "nct_id": "NCT00513019"},
        {"drug": "CETIRIZINE", "condition": "colorectal neoplasms", "nct_id": "NCT01722162"},
        {"drug": "ONDANSETRON", "condition": "schizophrenia", "nct_id": "NCT02614586"},
        {"drug": "EZETIMIBE", "condition": "chronic hepatitis d", "nct_id": "NCT03099278"},
        {"drug": "MYCOPHENOLATE", "condition": "refractory marginal zone lymphoma", "nct_id": "NCT04195633"},
        {"drug": "METHOTREXATE", "condition": "persistent pregnancy of unknown location", "nct_id": "NCT01800162"},
        {"drug": "CYSTEINE", "condition": "parkinson's disease (pd)", "nct_id": "NCT07093944"}
    ]
    
    # 2. 21 True Passes from V4 fresh audit
    v4_passes = [
        {"drug": "EMPAGLIFLOZIN", "condition": "acute heart failure", "nct_id": "NCT03554200"},
        {"drug": "ESOMEPRAZOLE", "condition": "helicobacter pylori infections", "nct_id": "NCT06811207"},
        {"drug": "METOCLOPRAMIDE", "condition": "headache, migraine", "nct_id": "NCT02972502"},
        {"drug": "DOXAZOSIN", "condition": "kidney stones", "nct_id": "NCT00713739"},
        {"drug": "EMPAGLIFLOZIN", "condition": "non-alcoholic fatty liver disease", "nct_id": "NCT02637973"},
        {"drug": "BARICITINIB", "condition": "major depressive disorder", "nct_id": "NCT07003997"},
        {"drug": "BAZEDOXIFENE", "condition": "schizophreniform disorders", "nct_id": "NCT04113993"},
        {"drug": "ACETAZOLAMIDE", "condition": "multiple sclerosis", "nct_id": "NCT02466074"},
        {"drug": "NALTREXONE", "condition": "breast cancer", "nct_id": "NCT00379197"},
        {"drug": "ATORVASTATIN", "condition": "tumor responses", "nct_id": "NCT03275376"},
        {"drug": "NALTREXONE", "condition": "melanoma", "nct_id": "NCT01650350"},
        {"drug": "ONDANSETRON", "condition": "alcohol abuse", "nct_id": "NCT00890149"},
        {"drug": "METRONIDAZOLE", "condition": "small intestinal bacterial overgrowth", "nct_id": "NCT02458781"},
        {"drug": "ACALABRUTINIB", "condition": "neuropathy;peripheral", "nct_id": "NCT05065554"},
        {"drug": "OXYTOCIN", "condition": "cervical ripening", "nct_id": "NCT03138252"},
        {"drug": "GLICLAZIDE", "condition": "psoriasis", "nct_id": "NCT01991197"},
        {"drug": "EZETIMIBE", "condition": "atherosclerosis", "nct_id": "NCT00651014"},
        {"drug": "TRYPTOPHAN", "condition": "alzheimer's disease", "nct_id": "NCT00202124"},
        {"drug": "CELECOXIB", "condition": "neuroblastoma", "nct_id": "NCT00357500"},
        {"drug": "ERDAFITINIB", "condition": "recurrent neuroblastoma", "nct_id": "NCT03210714"},
        {"drug": "DICLOFENAC", "condition": "diabetic macular edema", "nct_id": "NCT03458923"}
    ]
    
    # 3. 2 Checkpoint Passes
    checkpoint_passes = [
        {"drug": "VADADUSTAT", "condition": "acute respiratory distress syndrome", "nct_id": "NCT04478071"},
        {"drug": "EMPAGLIFLOZIN", "condition": "major depressive disorder", "nct_id": "NCT05757791"}
    ]
    
    # 4. Phase 1
    phase1_df = pd.read_csv('training_labels.csv')
    phase1_positives = phase1_df[phase1_df['label'] == 1]
    phase1_list = []
    for _, row in phase1_positives.iterrows():
        phase1_list.append({
            "drug": row['drug_name'],
            "condition": row['disease'],
            "nct_id": "Phase_1_Benchmark"
        })
    
    # Combine all
    all_raw = post_fix + v4_passes + checkpoint_passes + phase1_list
    
    # Deduplicate
    final_positives = []
    seen = set()
    for p in all_raw:
        pair = (p['drug'].upper(), p['condition'].lower())
        if pair not in seen:
            seen.add(pair)
            final_positives.append(p)
            
    with open('definitive_verified_positives.json', 'w') as f:
        json.dump(final_positives, f, indent=2)
        
    print(f"Total Combined: {len(all_raw)}")
    print(f"Total Deduplicated: {len(final_positives)}")
    
    # Save a text version for the user
    with open('definitive_list.txt', 'w') as f:
        f.write("--- POST-FIX AUDIT (9 CASES) ---\n")
        for x in post_fix:
            f.write(f"{x['drug']} -> {x['condition']} ({x['nct_id']})\n")
            
        f.write("\n--- V4 FRESH AUDIT (21 CASES) ---\n")
        for x in v4_passes:
            f.write(f"{x['drug']} -> {x['condition']} ({x['nct_id']})\n")
            
        f.write("\n--- CHECKPOINT (2 CASES) ---\n")
        for x in checkpoint_passes:
            f.write(f"{x['drug']} -> {x['condition']} ({x['nct_id']})\n")
            
        f.write("\n--- PHASE 1 BENCHMARKS ---\n")
        for x in phase1_list:
            f.write(f"{x['drug']} -> {x['condition']}\n")

if __name__ == "__main__":
    build()
