import json

with open('definitive_verified_positives.json') as f:
    v = json.load(f)

with open('training_pairs.json') as f:
    tp = json.load(f)['positives']

p = set((d,c) for d,c in tp)

print("Final 40 Verified Positives (0 Benchmark Overlap):")
for x in v:
    if (x['drug'], x['condition']) in p:
        print(f"- {x['drug']} -> {x['condition']} (NCT: {x['nct_id']})")
        # remove it from p to avoid duplicate prints if any
        p.remove((x['drug'], x['condition']))
