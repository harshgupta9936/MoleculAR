import json
import sqlite3
import random

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

def main():
    drugs = load_data()
    
    with open('definitive_verified_positives.json', 'r') as f:
        approved = json.load(f)
        
    benchmarks = {'SILDENAFIL', 'MINOXIDIL', 'METFORMIN', 'THALIDOMIDE', 'DULOXETINE', 'GABAPENTIN', 'BUPROPION', 'RALOXIFENE'}
    
    positives_raw = []
    for c in approved:
        dname = c['drug']
        cname = c['condition']
        if dname not in benchmarks and dname in drugs:
            positives_raw.append((dname, cname))
            
    # De-duplicate
    positives = list(set(positives_raw))
    print(f"Total valid unique positive pairs: {len(positives)}")
    
    # We will generate ~1.5x negatives
    num_negatives = int(len(positives) * 1.5)
    
    print("Generating negative pairs...")
    negatives = set()
    all_conditions = list(set([p[1] for p in positives]))
    drug_list = [d for d in drugs.keys() if d not in benchmarks and drugs[d]['indication'] != 'other' and drugs[d]['indication']]
    
    # For target overlap logic
    cond_to_tested_drugs = {}
    for (d, c) in positives:
        if c not in cond_to_tested_drugs: cond_to_tested_drugs[c] = set()
        cond_to_tested_drugs[c].add(d)
        
    attempts = 0
    while len(negatives) < num_negatives and attempts < 100000:
        attempts += 1
        dname = random.choice(drug_list)
        cname = random.choice(all_conditions)
        if (dname, cname) in positives: continue
        
        tested_drugs_for_c = cond_to_tested_drugs[cname]
        if dname in tested_drugs_for_c: continue
        
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
    
    with open('training_pairs.json', 'w') as f:
        json.dump({'positives': positives, 'negatives': negatives}, f, indent=2)
        
    overlap = [d for p in positives+negatives for d in p if d in benchmarks]
    print(f"Benchmark overlap: {len(overlap)}")

if __name__ == "__main__":
    main()
