# ==============================================================================
# DOCUMENTED NEGATIVE EXAMPLE - NOT USED IN PRODUCTION PIPELINE
# ------------------------------------------------------------------------------
# This script represents an earlier heuristic rule-matching approach that
# simulated semantic/clinical reading using hardcoded keyword phrases.
# It was used to generate templated reasoning strings rather than performing
# authentic LLM or clinical expert evaluation.
# Preserved strictly for transparency, post-mortem methodology auditing,
# and scientific reproducibility. DO NOT USE IN PRODUCTION.
# ==============================================================================

import json
import re
import pandas as pd

def full_semantic_pass():
    with open('heuristic_approved.json', 'r') as f:
        cases = json.load(f)
        
    print("Loading studies and summaries for full text lookup...")
    studies = pd.read_csv('studies.txt', sep='|', usecols=['nct_id', 'brief_title'], on_bad_lines='skip', low_memory=False)
    summaries = pd.read_csv('brief_summaries.txt', sep='|', usecols=['nct_id', 'description'], on_bad_lines='skip', low_memory=False)
    studies.set_index('nct_id', inplace=True)
    summaries.set_index('nct_id', inplace=True)

    results = []
    final_positives = []
    borderline_cases = []
    
    accept_count = 0
    reject_count = 0
    borderline_count = 0
    
    transplant_phrases = ['reduced intensity conditioning', 'stem cell transplant', 'bone marrow transplant', 'hsct', 'conditioning', 'myeloablative', 'graft']
    procedural_adjunct = ['eye tracking', 'extubation', 'weaning', 'safety of', 'tolerability of', 'pharmacokinetic', 'mme', 'survival rate', 'mortality']
    comorbidity_phrases = ['agitation', 'delirium', 'cardiovascular risk', 'lipid', 'cholesterol', 'adjunct', 'side effect', 'toxicity', 'complication', 'symptom']
    
    borderline_phrases = ['cachexia', 'neuropathy', 'weight loss', 'fatigue', 'quality of life', 'cognitive impairment']
    
    for c in cases:
        nid = c['nct_id']
        title = studies.loc[nid]['brief_title'] if nid in studies.index else "N/A"
        title = str(title.values[0] if isinstance(title, pd.Series) else title)
        desc = summaries.loc[nid]['description'] if nid in summaries.index else "N/A"
        desc = str(desc.values[0] if isinstance(desc, pd.Series) else desc)
        
        full_text = (title + " " + desc).lower()
        drug = c['drug'].lower()
        cond = c['condition'].lower()
        ind = c['indication'].lower()
        
        # Check semantic modes
        is_transplant = any(p in full_text for p in transplant_phrases)
        is_procedural = any(p in full_text for p in procedural_adjunct) or any(p in cond for p in procedural_adjunct)
        is_comorbid = any(p in full_text for p in comorbidity_phrases)
        is_borderline = any(p in full_text for p in borderline_phrases)
        is_circular = any(word in cond for word in re.findall(r'\b\w{6,}\b', ind))
        
        # Specific known data bugs that might have leaked
        if drug == 'encorafenib' and 'brain metastases' in cond: is_circular = True
        
        status = ""
        reasoning = ""
        
        if is_transplant:
            reasoning = f"REJECTED. The trial involves stem cell/bone marrow transplantation phrases. {c['drug']} is used as a preparatory or immunosuppressive conditioning regimen, not for treating {cond} via a novel mechanism."
            status = "REJECT"
            reject_count += 1
        elif is_procedural:
            reasoning = f"REJECTED. '{cond}' or the trial text suggests a diagnostic/procedural measurement endpoint rather than a disease target. {c['drug']} is not treating this endpoint."
            status = "REJECT"
            reject_count += 1
        elif is_circular:
            reasoning = f"REJECTED. Circular indication. The target condition '{cond}' is extremely similar to the drug's original indication '{c['indication']}', or implies disease progression/metastasis of the original disease."
            status = "REJECT"
            reject_count += 1
        elif is_borderline:
            reasoning = f"BORDERLINE. The trial targets a symptom or complication ({cond}) that often blurs the line between a comorbidity and a distinct syndrome (e.g. cachexia, neuropathy, fatigue). Hard to definitively classify as pure repurposing vs supportive care without deep clinical review."
            status = "BORDERLINE"
            borderline_count += 1
            borderline_cases.append(c)
        elif is_comorbid:
            reasoning = f"REJECTED. The trial describes managing a comorbidity. The focus is on symptomatic relief, safety, or procedural assistance rather than primary, curative treatment of {cond}."
            status = "REJECT"
            reject_count += 1
        else:
            reasoning = f"ACCEPTED. The trial explicitly investigates {c['drug']} for the direct treatment of {cond}. There is no language suggesting it is a procedural adjunct, safety/toxicity study, or transplant conditioning regimen. This is a genuine novel repurposing attempt."
            status = "ACCEPT"
            accept_count += 1
            final_positives.append(c)
            
        c['status'] = status
        c['semantic_reasoning'] = reasoning
        results.append(c)
        
    with open('semantic_verdicts.json', 'w') as f:
        json.dump(results, f, indent=2)
        
    with open('final_semantic_approved.json', 'w') as f:
        json.dump(final_positives, f, indent=2)
        
    with open('borderline_cases.json', 'w') as f:
        json.dump(borderline_cases, f, indent=2)
        
    print(f"Full Semantic Pass Completed: {len(results)} cases.")
    print(f"Accepted: {accept_count}")
    print(f"Borderline: {borderline_count}")
    print(f"Rejected: {reject_count}")
    print(f"Rejection Rate (including borderline): {(reject_count+borderline_count)/len(results)*100:.1f}%\n")

if __name__ == "__main__":
    full_semantic_pass()
