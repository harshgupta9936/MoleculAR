"""
experiments/step11_p4_protocol.py

Requirement N8:
Stricter Cross-Class Protocol P4 Definition & Query Count Audit:
- Relevant drugs share an indication in the primary gold standard but have NO ATC level-1 letter in common with query q:
    ATC_1(q) ∩ ATC_1(c) == ∅
- Candidates restricted to pool drugs with at least one ATC code.
- Duplicate-neutral: candidate pool excludes same-drug duplicates of q, and relevant set excludes duplicates of q.
- Computes exact query count, relevant pair count, pool size, analytical random expectation.
- Writes experiments/PROTOCOL_ADDENDUM.md.
- DOES NOT run M1 to M8 yet (waits for approval of the text).
"""

import sys
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from experiments.eval_engine import (
    load_database, load_gold_pairs, load_duplicates, rand_hit_at_k,
    DB_A_PATH, GOLD_PRIMARY_PATH, DATA_DIR, RESULTS_DIR
)

def load_atc_level1(pool_cids: set):
    rejected_sids = set()
    with open(DATA_DIR / "rejected_rule_b.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rejected_sids.add(int(r["dc_struct_id"]))

    import sqlite3
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

    cid_atc_l1 = defaultdict(set)
    cid_atc_full = defaultdict(set)
    with open(DATA_DIR / "drugcentral_raw_atc.csv", "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            sid = int(r["struct_id"])
            if sid in struct_to_cid:
                cid = struct_to_cid[sid]
                if cid in pool_cids:
                    code = r["atc_code"].strip() if "atc_code" in r and r["atc_code"] else ""
                    if code:
                        cid_atc_full[cid].add(code)
                        cid_atc_l1[cid].add(code[0].upper())

    return cid_atc_l1, cid_atc_full

def main():
    print("=" * 80)
    print("STEP 11: PROTOCOL P4 AUDIT & PROTOCOL_ADDENDUM.MD GENERATION")
    print("=" * 80)
    
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    drug_inds, _ = load_gold_pairs(GOLD_PRIMARY_PATH, set(pool_cids))
    duplicates_of, _ = load_duplicates(pool_cids)
    
    cid_atc_l1, cid_atc_full = load_atc_level1(set(pool_cids))
    pool_with_atc = sorted([cid for cid in pool_cids if len(cid_atc_l1[cid]) > 0])
    print(f"Total pool drugs: {len(pool_cids)}")
    print(f"Pool drugs with >= 1 ATC code: {len(pool_with_atc)}")
    
    cids_with_inds = sorted([cid for cid in pool_cids if len(drug_inds[cid]) > 0])
    cid_ind_sets = {cid: drug_inds[cid] for cid in cids_with_inds}
    
    # Audit P4 queries
    queries_p4 = []
    rel_pair_count = 0
    rand_hits = []
    
    for q_cid in cids_with_inds:
        if q_cid not in cid_atc_l1 or len(cid_atc_l1[q_cid]) == 0:
            continue
            
        q_l1 = cid_atc_l1[q_cid]
        q_inds = cid_ind_sets[q_cid]
        dups_q = duplicates_of.get(q_cid, set())
        
        # Base indication relevant in pool
        rel_P1 = [c_cid for c_cid in cids_with_inds if c_cid != q_cid and not q_inds.isdisjoint(cid_ind_sets[c_cid])]
        
        # P4 relevant: c in pool_with_atc, no common ATC level-1 letter, not in dups_q
        rel_p4 = [
            c for c in rel_P1
            if c in cid_atc_l1 and len(cid_atc_l1[c]) > 0
            and q_l1.isdisjoint(cid_atc_l1[c])
            and c not in dups_q
        ]
        
        if not rel_p4:
            continue
            
        # Candidates for q: all pool drugs with ATC, excluding q and dups_q
        cand_p4 = [c for c in pool_with_atc if c != q_cid and c not in dups_q]
        
        M_q = len(cand_p4)
        R_q = len(rel_p4)
        rel_pair_count += R_q
        rand_h10 = rand_hit_at_k(M_q, R_q, 10)
        rand_hits.append(rand_h10)
        
        queries_p4.append({
            "q_cid": q_cid,
            "q_name": drugs[q_cid]["name"],
            "q_atc_l1": sorted(list(q_l1)),
            "M": M_q,
            "R": R_q,
            "rand_hit10": rand_h10,
            "rel_cids": rel_p4
        })
        
    n_queries = len(queries_p4)
    analytical_rand_hit10 = float(np.mean(rand_hits)) if rand_hits else 0.0
    mean_M = float(np.mean([q['M'] for q in queries_p4]))
    mean_R = float(np.mean([q['R'] for q in queries_p4]))
    min_M = min(q['M'] for q in queries_p4)
    max_M = max(q['M'] for q in queries_p4)
    min_R = min(q['R'] for q in queries_p4)
    max_R = max(q['R'] for q in queries_p4)
    
    print(f"P4 Queries: {n_queries}")
    print(f"P4 Relevant Pairs: {rel_pair_count}")
    print(f"P4 Mean Candidate Pool Size (M): {mean_M:.1f}")
    print(f"P4 Mean Relevant Count (R): {mean_R:.2f}")
    print(f"P4 Analytical Random Hit@10: {analytical_rand_hit10:.6f} ({analytical_rand_hit10*100:.2f}%)")
    
    # Save query census to experiments/results/p4_query_census.csv
    census_csv = RESULTS_DIR / "p4_query_census.csv"
    with open(census_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["protocol", "P4"])
        w.writerow(["n_queries", n_queries])
        w.writerow(["n_relevant_pairs", rel_pair_count])
        w.writerow(["pool_drugs_with_atc", len(pool_with_atc)])
        w.writerow(["mean_M", f"{mean_M:.2f}"])
        w.writerow(["mean_R", f"{mean_R:.4f}"])
        w.writerow(["analytical_random_hit10", f"{analytical_rand_hit10:.6f}"])
        
    # Write PROTOCOL_ADDENDUM.md
    addendum_path = ROOT / "experiments" / "PROTOCOL_ADDENDUM.md"
    lines = [
        "# PROTOCOL ADDENDUM: Protocol P4 Specification",
        "",
        "**Date:** 2026-11-14",
        "**Status:** DRAFT - PENDING APPROVAL (Do not execute evaluation runs until explicitly approved)",
        f"**Parent Document:** [`experiments/PROTOCOL.md`](file:///{ROOT.as_posix()}/experiments/PROTOCOL.md) (Frozen - unmodified)",
        "",
        "---",
        "",
        "## 1. Motivation & Purpose",
        "",
        "Protocol P3 enforces cross-class drug repurposing by requiring that candidate and query drugs share no ATC Level-3 (pharmacological subgroup) code. However, drugs within the same anatomical or organ system (ATC Level 1) may still share broad physiological or target class overlap.",
        "",
        "**Protocol P4** establishes a stricter anatomical/organ-level cross-class barrier:",
        "A candidate drug is considered relevant to query drug $q$ if and only if they share a validated therapeutic indication in the gold standard but have **NO ATC Level-1 anatomical main group letter in common**.",
        "",
        "This evaluates true physiological cross-domain drug repurposing (e.g., repurposing a dermatological agent for nervous system disease, or an antineoplastic agent for cardiovascular disease).",
        "",
        "---",
        "",
        "## 2. Formal Protocol Specification",
        "",
        "### 2.1 Eligibility and Pools",
        "- **Drug Universe:** Primary pool $P$ consisting of $2,192$ non-trivial drugs from DB-A (excluding the 10 inorganic compounds defined in PROTOCOL.md Section 1.3).",
        "- **ATC Eligibility:** Drugs must have at least one valid WHO ATC code mapped via DrugCentral (`drugcentral_raw_atc.csv`), excluding rejected Rule (b) structures.",
        f"- Pool drugs with $\\ge 1$ ATC code: **{len(pool_with_atc)}** drugs.",
        "- **Candidate Pool for Query $q$ ($C_q$):**",
        "- Restricted strictly to pool drugs having at least one ATC code ($c \\in P_{\\text{ATC}}$).",
        "- Excludes the query drug itself ($c \\ne q$).",
        "- **Duplicate-Neutral:** Excludes any same-drug duplicate of $q$ (matching generic base name, matching first-block InChIKey, or matching desalted canonical SMILES).",
        "- Effective candidate pool size: $M_q = |P_{\\text{ATC}}| - 1 - |\\text{Duplicates}(q) \\cap P_{\\text{ATC}}|$.",
        "",
        "### 2.2 Relevance Definition",
        "A candidate drug $c \\in C_q$ is **relevant** to query $q$ ($c \\in \\text{Rel}_4(q)$) if and only if:",
        "1. $c$ and $q$ share at least one therapeutic indication in the primary gold standard (`gold_pairs_primary.csv`):",
        "   $$\\text{Ind}(q) \\cap \\text{Ind}(c) \\ne \\emptyset$$",
        "2. $c$ and $q$ have **no ATC Level-1 letter in common**:",
        "   $$\\text{ATC}_1(q) \\cap \\text{ATC}_1(c) = \\emptyset$$",
        "   where $\\text{ATC}_1(d) = \\{ \\text{code}[0] : \\text{code} \\in \\text{ATC}(d) \\}$.",
        "3. $c$ is not a same-drug duplicate of $q$ ($c \\notin \\text{Duplicates}(q)$).",
        "",
        "### 2.3 Query Set",
        "A pool drug $q \\in P_{\\text{ATC}}$ is an eligible **query** if and only if it has at least one relevant candidate under the condition above:",
        "$$|\\text{Rel}_4(q)| \\ge 1$$",
        "",
        "---",
        "",
        "## 3. Dataset & Protocol Parameters (Audit Counts)",
        "",
        f"- **Number of Qualifying Queries ($N$):** **{n_queries}**",
        f"- **Total Relevant Pairs:** **{rel_pair_count}**",
        f"- **Mean Candidate Pool Size ($M$):** **{mean_M:.2f}** (Range: {min_M} – {max_M})",
        f"- **Mean Relevant Drugs per Query ($R$):** **{mean_R:.2f}** (Range: {min_R} – {max_R})",
        f"- **Analytical Random Hit@10 Expectation:** **{analytical_rand_hit10:.6f}** ({analytical_rand_hit10*100:.2f}%)",
        "- **Ties & Precision:**",
        "- Rank tie-handling: Exact mid-rank average ranking: $r(c) = 1 + N_{>} + 0.5 \\times N_{=}$",
        "- Precision: IEEE-754 float64.",
        "- Bootstrap: 1,000 resamples (seed 42) for 95% confidence intervals.",
        ""
    ]
    with open(addendum_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
        
    print(f"PROTOCOL_ADDENDUM.md written to {addendum_path}")

if __name__ == "__main__":
    main()
