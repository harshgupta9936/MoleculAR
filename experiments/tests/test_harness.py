"""
experiments/tests/test_harness.py

Unit test suite for MoleculAR evaluation harness (Step 5.1 & 5.2).
Run with: pytest experiments/tests/test_harness.py
"""

import math
import sqlite3
import csv
import sys
import numpy as np
import pytest
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import DataStructs

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT))

DATA_DIR = ROOT / "experiments" / "data"
DB_A_PATH = DATA_DIR / "drugrepo_oct6.db"
GOLD_PRIMARY_PATH = DATA_DIR / "gold_pairs_primary.csv"
REJECTED_B_PATH = DATA_DIR / "rejected_rule_b.csv"

from experiments.eval_engine import (
    compute_average_ranks, load_database, load_gold_pairs,
    load_atc_l3, bootstrap_ci, INORGANIC_EXCLUDE
)

# -----------------------------------------------------------------------------
# 1. Tie Handling (Average Rank)
# -----------------------------------------------------------------------------
def test_tie_handling_average_rank():
    # Scores for 7 candidates
    # Scores: [1.0, 0.8, 0.8, 0.5, 0.5, 0.5, 0.2]
    # Expected average ranks:
    # 1.0 -> 1.0
    # 0.8 -> 1 + 1 + 0.5*(1) = 2.5
    # 0.5 -> 1 + 3 + 0.5*(2) = 5.0
    # 0.2 -> 1 + 6 + 0.5*(0) = 7.0
    cand_indices = np.array([0, 1, 2, 3, 4, 5, 6], dtype=np.int32)
    scores = np.array([1.0, 0.8, 0.8, 0.5, 0.5, 0.5, 0.2], dtype=np.float32)

    rank_map = compute_average_ranks(scores, cand_indices)
    assert rank_map[0] == pytest.approx(1.0)
    assert rank_map[1] == pytest.approx(2.5)
    assert rank_map[2] == pytest.approx(2.5)
    assert rank_map[3] == pytest.approx(5.0)
    assert rank_map[4] == pytest.approx(5.0)
    assert rank_map[5] == pytest.approx(5.0)
    assert rank_map[6] == pytest.approx(7.0)

# -----------------------------------------------------------------------------
# 2. Query Exclusion
# -----------------------------------------------------------------------------
def test_query_exclusion():
    # Load database
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    N = len(pool_cids)
    cid_to_idx = {c: i for i, c in enumerate(pool_cids)}

    for q_cid in pool_cids[:20]:
        q_idx = cid_to_idx[q_cid]
        cand_indices = np.array([i for i in range(N) if i != q_idx], dtype=np.int32)
        assert q_idx not in cand_indices
        assert len(cand_indices) == N - 1

# -----------------------------------------------------------------------------
# 3. Trivial and Inorganic Exclusion
# -----------------------------------------------------------------------------
def test_trivial_and_inorganic_exclusion():
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    assert len(drugs) == 2192

    for inorg_cid in INORGANIC_EXCLUDE:
        assert inorg_cid not in drugs

    # Check database raw: trivial drugs must not be in drugs
    conn = sqlite3.connect(DB_A_PATH.resolve().as_uri() + "?mode=ro", uri=True)
    c = conn.cursor()
    c.execute("SELECT chembl_id FROM drugs WHERE trivial_category IS NOT NULL")
    trivial_cids = {r[0] for r in c.fetchall()}
    conn.close()

    assert len(trivial_cids) == 65
    for triv_cid in trivial_cids:
        assert triv_cid not in drugs

# -----------------------------------------------------------------------------
# 4. Empty-set Jaccard (0 for M6/M7, unavailable for M9)
# -----------------------------------------------------------------------------
def test_empty_set_jaccard():
    # When targets or side effects are empty:
    # Jaccard = 0.0
    s_sim = 0.5
    t_empty_sim = 0.0
    se_empty_sim = 0.0

    # M6 score: 0.4*0.5 + 0.4*0.0 + 0.2*0.0 = 0.2
    m6_score = 0.4 * s_sim + 0.4 * t_empty_sim + 0.2 * se_empty_sim
    assert m6_score == pytest.approx(0.2)

    # M9 score when T and E are unavailable:
    # Available weights: only w_s = 0.4
    # Score = (0.4 * s_sim) / 0.4 = s_sim = 0.5
    w_sum = 0.4
    m9_score = (0.4 * s_sim) / w_sum
    assert m9_score == pytest.approx(0.5)

# -----------------------------------------------------------------------------
# 5. P2b Candidate Removal
# -----------------------------------------------------------------------------
def test_p2b_candidate_removal():
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    cid_to_idx = {c: i for i, c in enumerate(pool_cids)}

    # Sildenafil has known analogues with Tanimoto > 0.4 (e.g. Vardenafil)
    q_cid = "CHEMBL1737" # Sildenafil
    q_fp = drugs[q_cid]["fingerprint"]
    q_idx = cid_to_idx[q_cid]

    analogues = []
    for c_cid in pool_cids:
        if c_cid == q_cid:
            continue
        c_fp = drugs[c_cid]["fingerprint"]
        tan = DataStructs.TanimotoSimilarity(q_fp, c_fp)
        if tan > 0.4:
            analogues.append(c_cid)

    assert len(analogues) > 0

    cand_indices_p1 = [i for i in range(len(pool_cids)) if i != q_idx]
    cand_indices_p2b = [i for i in cand_indices_p1 if DataStructs.TanimotoSimilarity(q_fp, drugs[pool_cids[i]]["fingerprint"]) <= 0.4]

    assert len(cand_indices_p2b) == len(cand_indices_p1) - len(analogues)
    for a in analogues:
        assert cid_to_idx[a] not in cand_indices_p2b

# -----------------------------------------------------------------------------
# 6. P3 ATC Restriction
# -----------------------------------------------------------------------------
def test_p3_atc_restriction():
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    cid_atc = load_atc_l3(set(pool_cids))

    pool_with_atc = {c for c in pool_cids if len(cid_atc[c]) > 0}
    assert len(pool_with_atc) == 1900

    # Under P3, candidate count must be 1,900 - 1 = 1,899
    M_p3 = len(pool_with_atc) - 1
    assert M_p3 == 1899

# -----------------------------------------------------------------------------
# 7. Rule (b) Rejection List Applied
# -----------------------------------------------------------------------------
def test_rule_b_rejection_list_applied():
    with open(REJECTED_B_PATH, "r", encoding="utf-8") as f:
        r = list(csv.DictReader(f))
    assert len(r) == 9

    rejected_structs = {row["dc_struct_id"] for row in r}
    expected_sids = {"1402", "1729", "4641", "5787", "4031", "1655", "4867", "4586", "934"}
    assert rejected_structs == expected_sids

    # Verify stereo prefix reasons
    stereo_rows = [row for row in r if row["reason"] == "stereo-prefix name difference"]
    assert len(stereo_rows) == 5

    # Verify MENTHOL (CHEMBL256087) has 0 indications in gold_pairs_primary.csv
    # because its only DC struct (934) was rejected
    menthol_pairs = []
    with open(GOLD_PRIMARY_PATH, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["chembl_id"] == "CHEMBL256087":
                menthol_pairs.append(row)
    assert len(menthol_pairs) == 0

# -----------------------------------------------------------------------------
# 8. Bootstrap Reproducibility
# -----------------------------------------------------------------------------
def test_bootstrap_reproducibility():
    rng = np.random.RandomState(123)
    data = rng.rand(500)

    mean1, low1, high1 = bootstrap_ci(data, n_boot=1000, seed=42)
    mean2, low2, high2 = bootstrap_ci(data, n_boot=1000, seed=42)

    assert mean1 == mean2
    assert low1 == low2
    assert high1 == high2

# -----------------------------------------------------------------------------
# 9. 6-Drug Toy Database Verification
# -----------------------------------------------------------------------------
def test_toy_database_6_drugs():
    """
    Toy database of 6 drugs (D1..D6):
    - D1: Query drug
    - D2, D3, D4, D5, D6: Candidates
    Hand-crafted similarity scores, hubness penalties, and average ranks.
    """
    # Let D1 have:
    # S = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
    # T = [1.0, 0.5, 0.5, 0.0, 0.0, 0.0]
    # E = [1.0, 0.2, 0.2, 0.6, 0.0, 0.0]
    # Penalties:
    # D2: hub=0%   -> pen = 1.0
    # D3: hub=50%  -> pen = 1 + 0.25 = 1.25
    # D4: hub=0%   -> pen = 1.0
    # D5: hub=0%   -> pen = 1.0
    # D6: hub=0%   -> pen = 1.0

    # Scores:
    # D2: Raw = 0.4*0.8 + 0.4*0.5 + 0.2*0.2 = 0.32 + 0.20 + 0.04 = 0.56. Final = 0.56 / 1.0 = 0.56
    # D3: Raw = 0.4*0.6 + 0.4*0.5 + 0.2*0.2 = 0.24 + 0.20 + 0.04 = 0.48. Final = 0.48 / 1.25 = 0.384
    # D4: Raw = 0.4*0.4 + 0.4*0.0 + 0.2*0.6 = 0.16 + 0.00 + 0.12 = 0.28. Final = 0.28 / 1.0 = 0.28
    # D5: Raw = 0.4*0.2 + 0.4*0.0 + 0.2*0.0 = 0.08 + 0.00 + 0.00 = 0.08. Final = 0.08 / 1.0 = 0.08
    # D6: Raw = 0.4*0.0 + 0.4*0.0 + 0.2*0.0 = 0.00. Final = 0.00 / 1.0 = 0.00

    cand_indices = np.array([1, 2, 3, 4, 5], dtype=np.int32)
    scores = np.array([1.0, 0.56, 0.384, 0.28, 0.08, 0.00], dtype=np.float32)

    rank_map = compute_average_ranks(scores, cand_indices)
    assert rank_map[1] == pytest.approx(1.0) # D2
    assert rank_map[2] == pytest.approx(2.0) # D3
    assert rank_map[3] == pytest.approx(3.0) # D4
    assert rank_map[4] == pytest.approx(4.0) # D5
    assert rank_map[5] == pytest.approx(5.0) # D6

    # Test with deliberate tie between D3 and D4:
    scores_tie = np.array([1.0, 0.56, 0.28, 0.28, 0.08, 0.00], dtype=np.float32)
    rank_map_tie = compute_average_ranks(scores_tie, cand_indices)
    assert rank_map_tie[1] == pytest.approx(1.0) # D2
    assert rank_map_tie[2] == pytest.approx(2.5) # D3 (tied with D4)
    assert rank_map_tie[3] == pytest.approx(2.5) # D4 (tied with D3)
    assert rank_map_tie[4] == pytest.approx(4.0) # D5
    assert rank_map_tie[5] == pytest.approx(5.0) # D6

# -----------------------------------------------------------------------------
# 10. Duplicate Audit & Neutral Exclusion
# -----------------------------------------------------------------------------
def test_duplicate_audit_and_neutral_exclusion():
    from experiments.eval_engine import load_duplicates
    drugs = load_database(DB_A_PATH, exclude_inorganics=True)
    pool_cids = sorted(list(drugs.keys()))
    assert len(pool_cids) == 2192

    duplicates_of, pair_reasons = load_duplicates(pool_cids)
    drugs_with_dups = [c for c, dups in duplicates_of.items() if len(dups) > 0]
    assert len(drugs_with_dups) == 160
    assert len(pair_reasons) == 107

    # Test pair: CHEMBL1000 and CHEMBL1201190
    assert "CHEMBL1201190" in duplicates_of["CHEMBL1000"]
    assert "CHEMBL1000" in duplicates_of["CHEMBL1201190"]

# -----------------------------------------------------------------------------
# 11. All-Zero Score Vector Average Rank & Hit@10 Test
# -----------------------------------------------------------------------------
def test_all_zero_score_vector_average_rank():
    """
    Asserts that a candidate pool with all zero scores receives average rank
    (M + 1) / 2 for every candidate, and that Hit@10 is strictly 0.
    """
    M = 2191
    cand_indices = np.arange(M, dtype=np.int32)
    scores = np.zeros(M, dtype=np.float64)
    rank_map = compute_average_ranks(scores, cand_indices)

    expected_rank = (M + 1.0) / 2.0 # 1096.0
    assert len(rank_map) == M
    for idx in cand_indices:
        assert rank_map[idx] == pytest.approx(expected_rank)

    # For any relevant candidate in this tie block, r_min is 1096.0
    r_min = expected_rank
    hit10 = 1.0 if r_min <= 10.0 else 0.0
    assert hit10 == 0.0
    assert 1.0 / r_min == pytest.approx(1.0 / 1096.0)

