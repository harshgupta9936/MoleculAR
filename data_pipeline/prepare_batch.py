import json
import pandas as pd

def main():
    with open('heuristic_approved.json', 'r') as f:
        approved = json.load(f)
        
    studies = pd.read_csv('studies.txt', sep='|', usecols=['nct_id', 'brief_title'], on_bad_lines='skip', low_memory=False)
    summaries = pd.read_csv('brief_summaries.txt', sep='|', usecols=['nct_id', 'description'], on_bad_lines='skip', low_memory=False)
    studies.set_index('nct_id', inplace=True)
    summaries.set_index('nct_id', inplace=True)
    
    batch = []
    for i, s in enumerate(approved[:250]):
        nid = s['nct_id']
        title = studies.loc[nid]['brief_title'] if nid in studies.index else "N/A"
        title = str(title.values[0] if isinstance(title, pd.Series) else title)
        desc = summaries.loc[nid]['description'] if nid in summaries.index else "N/A"
        desc = str(desc.values[0] if isinstance(desc, pd.Series) else desc)[:200]
        
        batch.append({
            'id': i+1,
            'drug': s['drug'],
            'ind': s['indication'],
            'cond': s['condition'],
            'title': title,
            'desc': desc
        })
        
    with open('batch_250.json', 'w') as f:
        json.dump(batch, f, indent=2)

if __name__ == "__main__":
    main()
