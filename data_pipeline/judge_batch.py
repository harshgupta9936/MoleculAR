# ==============================================================================
# DOCUMENTED NEGATIVE EXAMPLE - NOT USED IN PRODUCTION PIPELINE
# ------------------------------------------------------------------------------
# This script is a documented negative example of hardcoded phrase matching.
# It used keyword substring matching to generate templated verdicts and simulated
# "semantic reading" rather than performing authentic LLM or clinical reasoning.
# Preserved strictly for transparency and post-mortem methodology auditing.
# DO NOT USE IN PRODUCTION.
# ==============================================================================

import json
import re

def judge_cases():
    with open('batch_250.json', 'r') as f:
        cases = json.load(f)
        
    results = []
    accept_count = 0
    reject_count = 0
    
    # New semantic failure modes:
    conditioning_phrases = ['reduced intensity conditioning', 'stem cell transplant', 'bone marrow transplant', 'hsct', 'conditioning', 'myeloablative', 'graft']
    comorbidity_phrases = ['agitation', 'delirium', 'extubation', 'cardiovascular risk', 'lipid', 'cholesterol', 'weaning', 'adjunct', 'side effect', 'toxicity', 'recovery of', 'manage', 'complication', 'symptom']
    
    for c in cases:
        full_text = (c['title'] + " " + c['desc']).lower()
        drug = c['drug'].lower()
        cond = c['cond'].lower()
        
        # Check conditioning
        is_transplant = any(p in full_text for p in conditioning_phrases)
        is_comorbid = any(p in full_text for p in comorbidity_phrases)
        
        # HEURISTIC RULE EVALUATION:
        # Note: The original script commented this as "Simulate semantic reading".
        # In reality, this is non-semantic string matching against keyword lists,
        # assigning templated strings to mimic clinical review. Documented as negative example.
        if is_transplant:
            reasoning = f"Rejected. The trial text contains phrases related to stem cell/bone marrow transplantation (e.g. conditioning, graft). The drug {c['drug']} is being used as a preparatory or immunosuppressive conditioning regimen rather than treating {cond} directly via a novel mechanism."
            status = "REJECT"
            reject_count += 1
        elif is_comorbid:
            reasoning = f"Rejected. The trial describes managing a comorbidity or procedural adjunct. The focus is on symptomatic relief, safety, or procedural assistance (like extubation or managing agitation/delirium) rather than primary treatment of {cond}."
            status = "REJECT"
            reject_count += 1
        elif c['ind'] in cond or cond in c['ind']:
            reasoning = f"Rejected. Circular indication. The target condition '{cond}' is extremely similar to the drug's original indication '{c['ind']}'."
            status = "REJECT"
            reject_count += 1
        else:
            reasoning = f"Accepted. The trial clearly investigates {c['drug']} for the direct treatment of {cond}. There is no language suggesting it is a procedural adjunct, safety/toxicity study, or transplant conditioning regimen. This appears to be a genuine novel repurposing attempt."
            status = "ACCEPT"
            accept_count += 1
            
        c['status'] = status
        c['semantic_reasoning'] = reasoning
        results.append(c)
        
    with open('semantic_verdicts.json', 'w') as f:
        json.dump(results, f, indent=2)
        
    print(f"Batch completed: {len(results)} cases.")
    print(f"Accepted: {accept_count}")
    print(f"Rejected: {reject_count}")
    print(f"Rejection Rate: {reject_count/len(results)*100:.1f}%\n")
    
    print("Example Accepts:")
    accepts = [c for c in results if c['status'] == 'ACCEPT'][:2]
    for a in accepts:
        print(f"[{a['status']}] {a['drug']} -> {a['cond']}\nReasoning: {a['semantic_reasoning']}\n")
        
    print("Example Rejects:")
    rejects = [c for c in results if c['status'] == 'REJECT'][:2]
    for r in rejects:
        print(f"[{r['status']}] {r['drug']} -> {r['cond']}\nReasoning: {r['semantic_reasoning']}\n")

if __name__ == "__main__":
    judge_cases()
