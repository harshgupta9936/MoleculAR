import sqlite3
import pandas as pd
from rdkit.Chem import DataStructs
from scoring_engine import repurposing_score_bulk

# DrugFeaturesCache: keys are chembl_id (primary key), with fallback lookup by drug_name
class DrugFeaturesCache(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name_to_cid = {}

    def __getitem__(self, key):
        if super().__contains__(key):
            return super().__getitem__(key)
        norm_key = str(key).strip().upper()
        if super().__contains__(norm_key):
            return super().__getitem__(norm_key)
        if norm_key in self.name_to_cid:
            cid = self.name_to_cid[norm_key]
            return super().__getitem__(cid)
        raise KeyError(key)

    def __contains__(self, key):
        if super().__contains__(key):
            return True
        norm_key = str(key).strip().upper()
        return super().__contains__(norm_key) or norm_key in self.name_to_cid

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

# Global Cache
drug_features = DrugFeaturesCache()
hub_percentile_ranks = {}
name_to_chembl = {}
df = None  # Kept for backward compatibility with server.py calling `df["drug_name"]`

def resolve_chembl_id(drug_identifier: str):
    """Resolve a chembl_id or drug name (base name or full name) to a valid chembl_id."""
    if not drug_identifier:
        return None
    q = str(drug_identifier).strip().upper()
    if dict.__contains__(drug_features, q):
        return q
    return drug_features.name_to_cid.get(q)

def init_db_cache():
    global drug_features, hub_percentile_ranks, name_to_chembl, df
    print("Loading SQLite database into memory cache (keyed by chembl_id)...")
    conn = sqlite3.connect("drugrepo.db")
    
    # Load basic DataFrame for fast lookups
    df = pd.read_sql("SELECT chembl_id, generic_base_name as drug_name, generic_base_name, indication, hubness_percentile FROM drugs", conn)
    
    cursor = conn.cursor()
    cursor.execute("SELECT chembl_id, generic_base_name, drug_name, fingerprint_blob, hubness_percentile, trivial_category FROM drugs")
    for row in cursor.fetchall():
        cid, base_name, full_name, fp_blob, hub_pct, trivial = row
        display_name = str(base_name).strip().upper() if base_name else str(full_name).strip().upper()
        fp = DataStructs.CreateFromBitString(fp_blob) if fp_blob else None
        
        # Primary keying by chembl_id
        drug_features[cid] = {
            "chembl_id": cid,
            "drug_name": display_name,
            "generic_base_name": display_name,
            "full_name": full_name,
            "fingerprint": fp,
            "targets": set(),
            "target_names": set(),
            "side_effects": set(),
            "hubness_percentile": hub_pct,
            "is_trivial": trivial is not None,
            "trivial_category": trivial
        }
        
        if hub_pct is not None:
            hub_percentile_ranks[cid] = hub_pct
            
        # Build name-to-chembl mapping
        if base_name:
            bn_upper = str(base_name).strip().upper()
            name_to_chembl[bn_upper] = cid
            drug_features.name_to_cid[bn_upper] = cid
        if full_name:
            fn_upper = str(full_name).strip().upper()
            name_to_chembl[fn_upper] = cid
            drug_features.name_to_cid[fn_upper] = cid
            
    # Load targets directly by chembl_id
    cursor.execute("SELECT chembl_id, uniprot_id, target_name FROM targets")
    for row in cursor.fetchall():
        cid, uid, tname = row
        if cid in drug_features:
            if uid:
                drug_features[cid]["targets"].add(uid)
            if tname:
                drug_features[cid]["target_names"].add(tname)
                
    # Load side effects directly by chembl_id
    cursor.execute("SELECT chembl_id, side_effect_name FROM side_effects")
    for row in cursor.fetchall():
        cid, se = row
        if cid in drug_features and se:
            drug_features[cid]["side_effects"].add(se)
            
    conn.close()
    print(f"Loaded {len(drug_features)} unique drugs (all keyed by chembl_id).")

# Initialize on import
init_db_cache()

def get_target_names_for_drug(drug_identifier: str):
    cid = resolve_chembl_id(drug_identifier)
    if cid and cid in drug_features:
        return drug_features[cid]["target_names"]
    return set()

def get_repurposing_candidates(drug_identifier: str, top_n=10, use_hub_correction=True):
    cid = resolve_chembl_id(drug_identifier)
    if not cid or cid not in drug_features:
        return []
        
    query_feats = drug_features[cid]
    query_name = query_feats["generic_base_name"]
    
    # Target drugs list: all drugs except query drug (and its salt duplicates) and trivial compounds
    target_drugs_list = []
    for other_cid, feats in drug_features.items():
        if other_cid != cid and feats.get("generic_base_name") != query_name and not feats.get("is_trivial", False):
            target_drugs_list.append(feats)
            
    # Bulk similarity score computation
    bulk_results = repurposing_score_bulk(query_feats, target_drugs_list, w_structure=0.4, w_target=0.4, w_sideeffect=0.2)
    
    scores = []
    for t_drug, score_info in zip(target_drugs_list, bulk_results):
        other_name = t_drug["generic_base_name"]
        raw_score = score_info["total_score"]
        
        hub_pct = t_drug.get("hubness_percentile")
        if use_hub_correction and hub_pct is not None:
            pct_rank = hub_pct / 100.0
            penalty_factor = 1.0 + (pct_rank ** 2)
            final_score = raw_score / penalty_factor
        else:
            final_score = raw_score
            
        scores.append((t_drug, other_name, final_score, score_info))
        
    scores.sort(key=lambda x: x[2], reverse=True)
    candidates = scores[:top_n]
    
    query_target_names = query_feats["target_names"]
    
    results = []
    for t_drug, candidate_name, score, score_info in candidates:
        candidate_target_names = t_drug["target_names"]
        shared_targets = query_target_names & candidate_target_names
        
        shared_ses = query_feats["side_effects"] & t_drug["side_effects"]
        top_shared_ses = sorted(shared_ses)[:5] if shared_ses else []
        
        ind = df[df["chembl_id"] == t_drug["chembl_id"]]["indication"].values
        ind_str = ind[0] if len(ind) > 0 and ind[0] else "Unknown"
        
        results.append({
            "candidate": candidate_name,
            "chembl_id": t_drug["chembl_id"],
            "score": score,
            "struct_sim": score_info["struct_sim"],
            "target_sim": score_info["target_sim"],
            "side_effect_sim": score_info["side_effect_sim"],
            "shared_targets": sorted(shared_targets) if shared_targets else [],
            "shared_side_effects": top_shared_ses,
            "indication": ind_str,
        })
        
    return results

def get_target_candidates(target_id: str, top_n=10, use_hub_correction=True):
    target_id = str(target_id).strip().upper()
    
    # Virtual drug feature vector for target search
    query_feats = {
        "chembl_id": f"TARGET_{target_id}",
        "drug_name": f"TARGET_{target_id}",
        "generic_base_name": f"TARGET_{target_id}",
        "fingerprint": None,
        "targets": {target_id},
        "target_names": {target_id},
        "side_effects": set()
    }
    
    target_drugs_list = list(drug_features.values())
    
    # Purely mechanistic target scoring
    bulk_results = repurposing_score_bulk(query_feats, target_drugs_list, w_structure=0.0, w_target=1.0, w_sideeffect=0.0)
    
    scores = []
    for t_drug, score_info in zip(target_drugs_list, bulk_results):
        if score_info["target_sim"] == 0:
            continue
            
        other_name = t_drug["generic_base_name"]
        raw_score = score_info["total_score"]
        
        # Safely handle NULL hubness so Mode B cannot raise TypeError
        hub_pct = t_drug.get("hubness_percentile")
        if use_hub_correction and hub_pct is not None:
            pct_rank = hub_pct / 100.0
            penalty_factor = 1.0 + (pct_rank ** 2)
            final_score = raw_score / penalty_factor
        else:
            final_score = raw_score
            
        scores.append((t_drug, other_name, final_score, score_info))
        
    if not scores:
        return []
        
    scores.sort(key=lambda x: x[2], reverse=True)
    candidates = scores[:top_n]
    
    results = []
    for t_drug, candidate_name, score, score_info in candidates:
        ind = df[df["chembl_id"] == t_drug["chembl_id"]]["indication"].values
        ind_str = ind[0] if len(ind) > 0 and ind[0] else "Unknown"
        
        results.append({
            "candidate": candidate_name,
            "chembl_id": t_drug["chembl_id"],
            "score": score,
            "struct_sim": score_info["struct_sim"],
            "target_sim": score_info["target_sim"],
            "side_effect_sim": score_info["side_effect_sim"],
            "shared_targets": [target_id],
            "shared_side_effects": [],
            "indication": ind_str,
        })
        
    return results
