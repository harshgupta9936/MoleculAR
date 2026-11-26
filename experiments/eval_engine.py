"""
experiments/eval_engine.py

Core Evaluation Engine for MoleculAR Benchmarks:
- Loads DB-A, DB-R, gold standards, ATC codes.
- Computes pairwise similarity matrices (S, T, E).
- Evaluates Methods M1 to M9 across Protocols P1, P2a, P2b, P3, P3-all.
- Average rank tie-handling: r(c) = 1 + N_> + 0.5 * N_=.
- Metrics: hit@10, hit@20, precision@10, MRR, percentile_rank, enrichment.
- 1,000 bootstrap resamples (seed 42) for 95% CIs with ratio inside each resample.
"""

import math
import sqlite3
import numpy as np
from collections import defaultdict, Counter
from pathlib import Path
import csv
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser

ROOT = Path(r"C:\Users\hg340\Documents\MoleculAR")
DATA_DIR = ROOT / "experiments" / "data"
RESULTS_DIR = ROOT / "experiments" / "results"
DB_A_PATH = DATA_DIR / "drugrepo_oct6.db"
DB_R_PATH = DATA_DIR / "drugrepo_oct6_recovered.db"
GOLD_PRIMARY_PATH = DATA_DIR / "gold_pairs_primary.csv"
EXCLUDED_IND_PATH = DATA_DIR / "excluded_indications.csv"
ATC_PATH = DATA_DIR / "drugcentral_raw_atc.csv"

INORGANIC_EXCLUDE = {
    "CHEMBL1200345", "CHEMBL1200458", "CHEMBL1200539", "CHEMBL1200691",
    "CHEMBL1200696", "CHEMBL1237066", "CHEMBL1356", "CHEMBL2068884",
    "CHEMBL2106119", "CHEMBL3833401",
}

def rand_hit_at_k(M: int, R: int, k: int = 10) -> float:
    if M < k:
        return 1.0 if R > 0 else 0.0
    if M - R < k:
        return 1.0
    return 1.0 - (math.comb(M - R, k) / math.comb(M, k))

def load_database(db_path: Path, exclude_inorganics: bool = True):
    conn = sqlite3.connect(db_path.resolve().as_uri() + "?mode=ro", uri=True)
    c = conn.cursor()
    c.execute("SELECT chembl_id, generic_base_name, trivial_category, fingerprint_blob, hubness_percentile FROM drugs")
    drugs = {}
    for cid, gbn, triv, fp_blob, hub in c.fetchall():
        if triv is not None:
            continue
        if exclude_inorganics and cid in INORGANIC_EXCLUDE:
            continue
        fp = DataStructs.CreateFromBitString(fp_blob) if fp_blob else None
        drugs[cid] = {
            "chembl_id": cid,
            "name": gbn,
            "fingerprint": fp,
            "hubness_percentile": hub if hub is not None else 0.0,
            "targets": set(),
            "side_effects": set(),
        }
    c.execute("SELECT chembl_id, uniprot_id FROM targets")
    for cid, uid in c.fetchall():
        if cid in drugs and uid:
            drugs[cid]["targets"].add(uid)
    c.execute("SELECT chembl_id, side_effect_name FROM side_effects")
    for cid, se in c.fetchall():
        if cid in drugs and se:
            drugs[cid]["side_effects"].add(se)
    conn.close()
    return drugs

def load_duplicates(pool_cids: list, db_path: Path = DB_A_PATH):
    conn = sqlite3.connect(db_path.resolve().as_uri() + "?mode=ro", uri=True)
    c = conn.cursor()
    c.execute("SELECT chembl_id, generic_base_name, smiles, trivial_category FROM drugs")
    pool_drugs = {}
    for cid, gbn, smiles, triv in c.fetchall():
        if triv is not None or cid in INORGANIC_EXCLUDE:
            continue
        pool_drugs[cid] = {
            "cid": cid,
            "name": gbn.strip() if gbn else "",
            "name_lower": gbn.strip().lower() if gbn else "",
            "smiles": smiles.strip() if smiles else "",
        }
    conn.close()

    lfc = LargestFragmentChooser()
    for cid, d in pool_drugs.items():
        mol = Chem.MolFromSmiles(d["smiles"])
        if mol:
            ik = Chem.MolToInchiKey(mol)
            d["ik_b1"] = ik.split("-")[0] if ik else ""
            mol_lfc = lfc.choose(mol)
            d["desalted_smiles"] = Chem.MolToSmiles(mol_lfc, canonical=True)
        else:
            d["ik_b1"] = ""
            d["desalted_smiles"] = ""

    duplicates_of = defaultdict(set)
    pair_reasons = {}

    for i in range(len(pool_cids)):
        c1 = pool_cids[i]
        d1 = pool_drugs.get(c1)
        if not d1: continue
        for j in range(i + 1, len(pool_cids)):
            c2 = pool_cids[j]
            d2 = pool_drugs.get(c2)
            if not d2: continue
            reasons = []
            if d1["name_lower"] and d1["name_lower"] == d2["name_lower"]:
                reasons.append("generic_base_name")
            if d1["ik_b1"] and d1["ik_b1"] == d2["ik_b1"]:
                reasons.append("first_block_inchikey")
            if d1["desalted_smiles"] and d1["desalted_smiles"] == d2["desalted_smiles"]:
                reasons.append("desalted_canonical_smiles")
            if reasons:
                duplicates_of[c1].add(c2)
                duplicates_of[c2].add(c1)
                pair_reasons[(c1, c2)] = reasons

    return duplicates_of, pair_reasons

def load_gold_pairs(gold_file: Path, pool_cids: set):
    drug_inds = defaultdict(set)
    ind_drugs = defaultdict(set)
    with open(gold_file, "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            cid = row["chembl_id"]
            ind = row["indication_id"]
            if cid in pool_cids:
                drug_inds[cid].add(ind)
                ind_drugs[ind].add(cid)
    return drug_inds, ind_drugs

def load_atc_l3(pool_cids: set):
    rejected_sids = set()
    with open(DATA_DIR / "rejected_rule_b.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rejected_sids.add(int(r["dc_struct_id"]))

    conn = sqlite3.connect(DB_A_PATH.resolve().as_uri() + "?mode=ro", uri=True)
    ik_full = {}
    ik_b1 = defaultdict(list)
    drugs_all = {}
    for cid, gbn, smiles, triv in conn.execute("SELECT chembl_id, generic_base_name, smiles, trivial_category FROM drugs").fetchall():
        drugs_all[cid] = {"triv": triv is not None}
        if smiles:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                ik = Chem.MolToInchiKey(mol)
                if ik:
                    ik_full[ik] = cid
                    ik_b1[ik.split("-")[0]].append(cid)
    cw_db, cw_pc = {}, {}
    for cid, db_id, pc_cid in conn.execute("SELECT chembl_id, drugbank_id, pubchem_cid FROM crosswalk").fetchall():
        if db_id and db_id.strip(): cw_db[db_id.strip()] = cid
        if pc_cid and str(pc_cid).strip(): cw_pc[str(pc_cid).strip()] = cid
    conn.close()

    import psycopg2
    conn_pg = psycopg2.connect(host="unmtid-dbs.net", port=5433, user="drugman", password="dosage", dbname="drugcentral", connect_timeout=30)
    cur_pg = conn_pg.cursor()
    cur_pg.execute("SELECT id, inchikey FROM structures ORDER BY id")
    all_structs = cur_pg.fetchall()
    cur_pg.execute("SELECT struct_id, id_type, identifier FROM identifier WHERE id_type IN ('DRUGBANK_ID', 'PUBCHEM_CID', 'ChEMBL_ID')")
    dc_idents = defaultdict(lambda: defaultdict(list))
    for sid, itype, ident in cur_pg.fetchall():
        if ident and ident.strip():
            dc_idents[sid][itype].append(ident.strip())
    conn_pg.close()

    struct_to_cid = {}
    for sid, ik in all_structs:
        if sid in rejected_sids:
            continue
        clean_ik = ik.strip() if ik else None
        b1 = clean_ik.split("-")[0] if clean_ik else None
        cid = None
        if clean_ik and clean_ik in ik_full:
            cid = ik_full[clean_ik]
        elif clean_ik and b1 in ik_b1:
            cid = ik_b1[b1][0]
        elif sid in dc_idents:
            ids = dc_idents[sid]
            for db_id in ids.get("DRUGBANK_ID", []):
                if db_id in cw_db: cid = cw_db[db_id]; break
            if not cid:
                for pc_id in ids.get("PUBCHEM_CID", []):
                    if pc_id in cw_pc: cid = cw_pc[pc_id]; break
            if not cid:
                for ch_id in ids.get("ChEMBL_ID", []):
                    if ch_id in drugs_all: cid = ch_id; break
        if cid:
            struct_to_cid[sid] = cid

    cid_atc_l3 = defaultdict(set)
    with open(ATC_PATH, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sid = int(r["struct_id"])
            if sid in struct_to_cid:
                cid = struct_to_cid[sid]
                if cid in pool_cids and r["l3_code"] and r["l3_code"].strip():
                    cid_atc_l3[cid].add(r["l3_code"].strip())

    return cid_atc_l3

def precompute_matrices(drugs: dict, pool_cids: list):
    N = len(pool_cids)
    cid_to_idx = {cid: i for i, cid in enumerate(pool_cids)}

    fps = [drugs[cid]["fingerprint"] for cid in pool_cids]
    S = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        fp_i = fps[i]
        if fp_i is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp_i, fps)
        S[i] = np.array(sims, dtype=np.float64)

    T = np.zeros((N, N), dtype=np.float64)
    has_target = np.array([len(drugs[cid]["targets"]) > 0 for cid in pool_cids], dtype=bool)
    for i in range(N):
        t_i = drugs[pool_cids[i]]["targets"]
        if not t_i:
            continue
        for j in range(i, N):
            t_j = drugs[pool_cids[j]]["targets"]
            if not t_j:
                continue
            u = len(t_i | t_j)
            val = len(t_i & t_j) / u if u > 0 else 0.0
            T[i, j] = val
            T[j, i] = val

    E = np.zeros((N, N), dtype=np.float64)
    has_se = np.array([len(drugs[cid]["side_effects"]) > 0 for cid in pool_cids], dtype=bool)
    for i in range(N):
        se_i = drugs[pool_cids[i]]["side_effects"]
        if not se_i:
            continue
        for j in range(i, N):
            se_j = drugs[pool_cids[j]]["side_effects"]
            if not se_j:
                continue
            u = len(se_i | se_j)
            val = len(se_i & se_j) / u if u > 0 else 0.0
            E[i, j] = val
            E[j, i] = val

    hub_pcts = np.array([drugs[cid]["hubness_percentile"] for cid in pool_cids], dtype=np.float64)
    penalties = 1.0 + (hub_pcts / 100.0) ** 2

    return {
        "S": S,
        "T": T,
        "E": E,
        "has_target": has_target,
        "has_se": has_se,
        "penalties": penalties,
        "hub_pcts": hub_pcts,
        "cid_to_idx": cid_to_idx,
    }

def compute_average_ranks(scores: np.ndarray, candidate_indices: np.ndarray):
    cand_scores = scores[candidate_indices]
    desc_unique, desc_counts = np.unique(cand_scores, return_counts=True)
    desc_order = np.argsort(-desc_unique)
    desc_unique = desc_unique[desc_order]
    desc_counts = desc_counts[desc_order]

    cum_counts = np.cumsum(desc_counts)
    n_gt = np.empty_like(desc_counts)
    n_gt[0] = 0
    if len(desc_counts) > 1:
        n_gt[1:] = cum_counts[:-1]

    avg_rank_per_unique = 1.0 + n_gt + 0.5 * (desc_counts - 1)
    idx_in_unique = np.searchsorted(-desc_unique, -cand_scores)
    cand_ranks = avg_rank_per_unique[idx_in_unique]

    return dict(zip(candidate_indices, cand_ranks))

def compute_m1_random(M: int, R: int, rng: np.random.RandomState, n_permutations: int = 100):
    hit10 = rand_hit_at_k(M, R, 10)
    hit20 = rand_hit_at_k(M, R, 20)
    prec10 = R / M

    min_ranks = []
    for _ in range(n_permutations):
        drawn = rng.choice(M, size=R, replace=False) + 1
        min_ranks.append(np.min(drawn))

    mrr = np.mean([1.0 / r for r in min_ranks])
    pct_rank = np.mean([r / M for r in min_ranks])

    return {
        "hit10": hit10,
        "hit20": hit20,
        "prec10": prec10,
        "mrr": mrr,
        "pct_rank": pct_rank,
        "rand_hit10": hit10,
    }

def bootstrap_ci(values: np.ndarray, n_boot: int = 1000, seed: int = 42):
    rng = np.random.RandomState(seed)
    n = len(values)
    boot_means = np.empty(n_boot, dtype=np.float64)
    indices = rng.randint(0, n, size=(n_boot, n))
    for b in range(n_boot):
        boot_means[b] = np.mean(values[indices[b]])
    ci_low = np.percentile(boot_means, 2.5)
    ci_high = np.percentile(boot_means, 97.5)
    return float(np.mean(values)), float(ci_low), float(ci_high)

def bootstrap_enrichment(hit10: np.ndarray, rand_hit10: np.ndarray, n_boot: int = 1000, seed: int = 42):
    rng = np.random.RandomState(seed)
    n = len(hit10)
    boot_enrich = np.empty(n_boot, dtype=np.float64)
    indices = rng.randint(0, n, size=(n_boot, n))
    for b in range(n_boot):
        idx = indices[b]
        m_hit = np.mean(hit10[idx])
        m_rand = np.mean(rand_hit10[idx])
        boot_enrich[b] = m_hit / m_rand if m_rand > 0 else 1.0
    mean_val = np.mean(hit10) / np.mean(rand_hit10) if np.mean(rand_hit10) > 0 else 1.0
    ci_low = np.percentile(boot_enrich, 2.5)
    ci_high = np.percentile(boot_enrich, 97.5)
    return float(mean_val), float(ci_low), float(ci_high)

def evaluate_protocol(
    protocol_name: str,
    pool_cids: list,
    matrices: dict,
    drug_inds: dict,
    cid_atc_l3: dict,
    seed: int = 42,
    theta: float = 0.4,
    duplicates_of: dict = None,
    k_hub: float = 2.0,
):
    """
    Evaluates Methods M1 to M9 for the given protocol.
    Returns per-method metrics and query details.
    """
    cid_to_idx = matrices["cid_to_idx"]
    S = matrices["S"]
    T = matrices["T"]
    E = matrices["E"]
    has_target = matrices["has_target"]
    has_se = matrices["has_se"]
    if k_hub == 2.0:
        penalties = matrices["penalties"]
    else:
        penalties = 1.0 + (matrices["hub_pcts"] / 100.0) ** k_hub
    N = len(pool_cids)

    if protocol_name in ("P1_dn", "P2a_dn_0.4", "P2b_dn_0.4", "P3_dn") and duplicates_of is None:
        duplicates_of, _ = load_duplicates(pool_cids)

    # Indication popularity scores for M8
    pop_scores = np.array([len(drug_inds[cid]) for cid in pool_cids], dtype=np.float64)

    # 1. Base P1 relevance
    # Precompute indication sets as frozensets
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}

    # Pool with ATC for P3
    pool_with_atc = {cid for cid in pool_cids if len(cid_atc_l3[cid]) > 0}
    atc_indices = np.array([cid_to_idx[cid] for cid in pool_with_atc], dtype=np.int32)

    queries = []
    # Identify valid queries and their relevant candidates
    for q_cid in cids_with_inds:
        q_idx = cid_to_idx[q_cid]
        q_inds = cid_ind_sets[q_cid]

        # P1 relevant
        rel_P1 = [c_cid for c_cid in cids_with_inds if c_cid != q_cid and not q_inds.isdisjoint(cid_ind_sets[c_cid])]
        if not rel_P1:
            continue

        if protocol_name == "P1":
            rel_cids = rel_P1
            cand_indices = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)

        elif protocol_name == "P1_dn":
            dups_q = duplicates_of.get(q_cid, set()) if duplicates_of else set()
            rel_cids = [c for c in rel_P1 if c not in dups_q]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in range(N) if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)

        elif protocol_name == "P2a_dn_0.4":
            dups_q = duplicates_of.get(q_cid, set()) if duplicates_of else set()
            rel_cids = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4 and c not in dups_q]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in range(N) if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)

        elif protocol_name.startswith("P2a"):
            # Analogue neutralization: filter relevant set, keep candidates full
            rel_cids = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= theta]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)

        elif protocol_name == "P2b_0.4":
            # Analogue removal: filter relevant set AND remove analogues from candidate pool
            rel_cids = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in range(N) if i != q_idx and S[q_idx, i] <= 0.4], dtype=np.int32)

        elif protocol_name == "P2b_dn_0.4":
            dups_q = duplicates_of.get(q_cid, set()) if duplicates_of else set()
            rel_cids = [c for c in rel_P1 if S[q_idx, cid_to_idx[c]] <= 0.4 and c not in dups_q]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in range(N) if i != q_idx and S[q_idx, i] <= 0.4 and pool_cids[i] not in dups_q], dtype=np.int32)

        elif protocol_name == "P3":
            # Primary P3: candidates restricted to pool drugs with >= 1 ATC code
            if q_cid not in pool_with_atc:
                continue
            q_atc = cid_atc_l3[q_cid]
            rel_cids = [c for c in rel_P1 if c in pool_with_atc and q_atc.isdisjoint(cid_atc_l3[c])]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in atc_indices if i != q_idx], dtype=np.int32)

        elif protocol_name == "P3_dn":
            if q_cid not in pool_with_atc:
                continue
            q_atc = cid_atc_l3[q_cid]
            dups_q = duplicates_of.get(q_cid, set()) if duplicates_of else set()
            rel_cids = [c for c in rel_P1 if c in pool_with_atc and q_atc.isdisjoint(cid_atc_l3[c]) and c not in dups_q]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in atc_indices if i != q_idx and pool_cids[i] not in dups_q], dtype=np.int32)

        elif protocol_name == "P3_all":
            # P3-all: relevant set same as P3, but all pool drugs are candidates
            if q_cid not in pool_with_atc:
                continue
            q_atc = cid_atc_l3[q_cid]
            rel_cids = [c for c in rel_P1 if c in pool_with_atc and q_atc.isdisjoint(cid_atc_l3[c])]
            if not rel_cids:
                continue
            cand_indices = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)

        else:
            raise ValueError(f"Unknown protocol {protocol_name}")

        rel_indices = np.array([cid_to_idx[c] for c in rel_cids], dtype=np.int32)
        queries.append({
            "q_cid": q_cid,
            "q_idx": q_idx,
            "rel_indices": rel_indices,
            "cand_indices": cand_indices,
            "M": len(cand_indices),
            "R": len(rel_indices),
        })

    n_queries = len(queries)
    if n_queries == 0:
        return {}, []

    # Container for query-level metrics per method
    method_names = ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"]
    query_metrics = {m: {
        "hit10": np.empty(n_queries, dtype=np.float64),
        "hit20": np.empty(n_queries, dtype=np.float64),
        "prec10": np.empty(n_queries, dtype=np.float64),
        "mrr": np.empty(n_queries, dtype=np.float64),
        "pct_rank": np.empty(n_queries, dtype=np.float64),
        "rand_hit10": np.empty(n_queries, dtype=np.float64),
    } for m in method_names}

    rng_m1 = np.random.RandomState(seed)

    for q_num, q in enumerate(queries):
        q_idx = q["q_idx"]
        rel_idx = q["rel_indices"]
        cand_idx = q["cand_indices"]
        M = q["M"]
        R = q["R"]

        # 1. M1 Random
        m1_res = compute_m1_random(M, R, rng_m1, n_permutations=100)
        for k in ["hit10", "hit20", "prec10", "mrr", "pct_rank", "rand_hit10"]:
            query_metrics["M1"][k][q_num] = m1_res[k]

        rand_hit10_q = m1_res["hit10"]

        # Precompute signal vectors for query q_idx
        s_row = S[q_idx]
        t_row = T[q_idx]
        e_row = E[q_idx]

        # M2 (S)
        scores_m2 = s_row

        # M3 (T)
        scores_m3 = t_row

        # M4 (E)
        scores_m4 = e_row

        # M5 (0.5 S + 0.5 T)
        scores_m5 = 0.5 * s_row + 0.5 * t_row

        # M6 (0.4 S + 0.4 T + 0.2 E)
        scores_m6 = 0.4 * s_row + 0.4 * t_row + 0.2 * e_row

        # M7 (M6 / penalty)
        scores_m7 = scores_m6 / penalties

        # M8 (Popularity)
        scores_m8 = pop_scores

        # M9 (Available signals reweighted / penalty)
        # S is always available (w=0.4)
        # T is available if has_target[q_idx] and has_target (w=0.4)
        # E is available if has_se[q_idx] and has_se (w=0.2)
        w_s = 0.4
        w_t = 0.4 if has_target[q_idx] else 0.0
        w_e = 0.2 if has_se[q_idx] else 0.0

        w_t_vec = w_t * has_target
        w_e_vec = w_e * has_se
        w_sum = w_s + w_t_vec + w_e_vec
        scores_m9_raw = (w_s * s_row + w_t_vec * t_row + w_e_vec * e_row) / w_sum
        scores_m9 = scores_m9_raw / penalties

        method_scores = {
            "M2": scores_m2,
            "M3": scores_m3,
            "M4": scores_m4,
            "M5": scores_m5,
            "M6": scores_m6,
            "M7": scores_m7,
            "M8": scores_m8,
            "M9": scores_m9,
        }

        for m_name, scores in method_scores.items():
            rank_map = compute_average_ranks(scores, cand_idx)
            rel_ranks = [rank_map[idx] for idx in rel_idx]
            r_min = np.min(rel_ranks)

            hit10 = 1.0 if r_min <= 10.0 else 0.0
            hit20 = 1.0 if r_min <= 20.0 else 0.0
            mrr = 1.0 / r_min
            prec10 = sum(1 for r in rel_ranks if r <= 10.0) / 10.0
            pct_rank = r_min / M

            qm = query_metrics[m_name]
            qm["hit10"][q_num] = hit10
            qm["hit20"][q_num] = hit20
            qm["prec10"][q_num] = prec10
            qm["mrr"][q_num] = mrr
            qm["pct_rank"][q_num] = pct_rank
            qm["rand_hit10"][q_num] = rand_hit10_q

    # Summarize with 1,000 bootstrap CIs
    summary_rows = []
    for m in method_names:
        qm = query_metrics[m]
        # hit@10
        m_hit10, low_hit10, high_hit10 = bootstrap_ci(qm["hit10"], n_boot=1000, seed=seed)
        summary_rows.append({"method": m, "n_queries": n_queries, "metric": "hit@10", "mean": m_hit10, "ci_low": low_hit10, "ci_high": high_hit10})

        # hit@20
        m_hit20, low_hit20, high_hit20 = bootstrap_ci(qm["hit20"], n_boot=1000, seed=seed)
        summary_rows.append({"method": m, "n_queries": n_queries, "metric": "hit@20", "mean": m_hit20, "ci_low": low_hit20, "ci_high": high_hit20})

        # precision@10
        m_prec10, low_prec10, high_prec10 = bootstrap_ci(qm["prec10"], n_boot=1000, seed=seed)
        summary_rows.append({"method": m, "n_queries": n_queries, "metric": "precision@10", "mean": m_prec10, "ci_low": low_prec10, "ci_high": high_prec10})

        # MRR
        m_mrr, low_mrr, high_mrr = bootstrap_ci(qm["mrr"], n_boot=1000, seed=seed)
        summary_rows.append({"method": m, "n_queries": n_queries, "metric": "mrr", "mean": m_mrr, "ci_low": low_mrr, "ci_high": high_mrr})

        # percentile rank
        m_pct, low_pct, high_pct = bootstrap_ci(qm["pct_rank"], n_boot=1000, seed=seed)
        summary_rows.append({"method": m, "n_queries": n_queries, "metric": "percentile_rank", "mean": m_pct, "ci_low": low_pct, "ci_high": high_pct})

        # enrichment
        m_enr, low_enr, high_enr = bootstrap_enrichment(qm["hit10"], qm["rand_hit10"], n_boot=1000, seed=seed)
        summary_rows.append({"method": m, "n_queries": n_queries, "metric": "enrichment", "mean": m_enr, "ci_low": low_enr, "ci_high": high_enr})

    return summary_rows, queries, query_metrics
