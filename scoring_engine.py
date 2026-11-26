from rdkit.Chem import DataStructs

def bulk_tanimoto(query_fp, target_fps_list):
    """Calculate Tanimoto similarity for a query against a list of fingerprints efficiently."""
    if query_fp is None or not target_fps_list:
        return [0.0] * len(target_fps_list)
    return DataStructs.BulkTanimotoSimilarity(query_fp, target_fps_list)

def jaccard_similarity(set1, set2):
    """Calculate Jaccard similarity for two sets (targets or side effects)."""
    set1, set2 = set(set1), set(set2)
    if not set1 or not set2:
        return 0.0
    return len(set1 & set2) / len(set1 | set2)

def repurposing_score_bulk(query_feats, target_drugs_list, w_structure=0.4, w_target=0.4, w_sideeffect=0.2):
    """
    Calculate combined similarity score in bulk.
    
    query_feats: { 'fingerprint': fp, 'targets': set(...), 'side_effects': set(...) }
    target_drugs_list: list of dicts with same keys
    
    Returns list of dicts with scores matching the input list order.
    """
    target_fps = [d.get("fingerprint") for d in target_drugs_list]
    struct_sims = bulk_tanimoto(query_feats.get("fingerprint"), target_fps)
    
    q_targs = query_feats.get("targets", set())
    q_ses = query_feats.get("side_effects", set())
    
    results = []
    for i, t_drug in enumerate(target_drugs_list):
        s_sim = struct_sims[i]
        t_sim = jaccard_similarity(q_targs, t_drug.get("targets", set()))
        se_sim = jaccard_similarity(q_ses, t_drug.get("side_effects", set()))
        
        total_score = (w_structure * s_sim) + (w_target * t_sim) + (w_sideeffect * se_sim)
        
        results.append({
            'total_score': total_score,
            'struct_sim': s_sim,
            'target_sim': t_sim,
            'side_effect_sim': se_sim
        })
        
    return results
