import json
import random

def extract():
    with open('semantic_verdicts.json', 'r') as f:
        cases = json.load(f)
        
    targets = {
        'ketamine': 'cancer',
        'oritavancin': 'substance use disorders',
        'olanzapine': 'germ cell tumor',
        'mirtazapine': "alzheimer's disease" # wait, the cond might be alzheimer disease
    }
    
    print("TARGETS:")
    for c in cases:
        drug = c['drug'].lower()
        cond = c['condition'].lower()
        if drug == 'ketamine' and 'cancer' in cond:
            print(f"[{c['status']}] {c['drug']} -> {c['condition']}\nReasoning: {c['semantic_reasoning']}\n")
        elif drug == 'oritavancin' and 'substance' in cond:
            print(f"[{c['status']}] {c['drug']} -> {c['condition']}\nReasoning: {c['semantic_reasoning']}\n")
        elif drug == 'olanzapine' and 'germ' in cond:
            print(f"[{c['status']}] {c['drug']} -> {c['condition']}\nReasoning: {c['semantic_reasoning']}\n")
        elif drug == 'mirtazapine' and 'alzheimer' in cond:
            print(f"[{c['status']}] {c['drug']} -> {c['condition']}\nReasoning: {c['semantic_reasoning']}\n")
            
    print("\n10 RANDOM ENTRIES (BACK HALF):")
    back_half = cases[len(cases)//2:]
    random.seed(42)
    sample = random.sample(back_half, 10)
    for c in sample:
        print(f"[{c['status']}] {c['drug']} -> {c['condition']}\nReasoning: {c['semantic_reasoning']}\n")

if __name__ == "__main__":
    extract()
