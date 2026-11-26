import pandas as pd
import json

with open('heuristic_approved_v8.json') as f:
    pos = json.load(f)

with open('generated_negatives.json') as f:
    neg = json.load(f)

data = []
for p in pos:
    data.append({
        'drug_name': p['drug'],
        'disease': p['condition'],
        'label': 1,
        'source': 'AACT Mining (V8)'
    })

for n in neg:
    data.append({
        'drug_name': n['drug'],
        'disease': n['condition'],
        'label': 0,
        'source': 'Generated Negative (V2)'
    })

df = pd.DataFrame(data)
df.to_csv('training_labels_v2.csv', index=False)
print(f"Saved {len(df)} pairs to training_labels_v2.csv ({len(pos)} pos, {len(neg)} neg)")
