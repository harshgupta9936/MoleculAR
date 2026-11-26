import os
import sys

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import sqlite3
import pandas as pd
from rdkit.Chem import DataStructs
from pipeline import get_repurposing_candidates, drug_features
from benchmark_cases import BENCHMARK_CASES_DICT

def test_chirality():
    print("\n--- Test 1: Chirality (Ofloxacin vs Levofloxacin) ---")
    if "OFLOXACIN" not in drug_features or "LEVOFLOXACIN" not in drug_features:
        print("Missing Ofloxacin or Levofloxacin in DB!")
        return False
        
    fp_o = drug_features["OFLOXACIN"]["fingerprint"]
    fp_l = drug_features["LEVOFLOXACIN"]["fingerprint"]
    
    sim = DataStructs.TanimotoSimilarity(fp_o, fp_l)
    print(f"Structural Similarity: {sim:.4f}")
    if sim < 1.0:
        print("PASS: Similarity dropped below 1.0, proving stereoisomer distinction.")
        return True
    else:
        print("FAIL: Similarity is still 1.0!")
        return False

def test_benchmark():
    print("\n--- Test 2: 8-Case Benchmark ---")
    cases = BENCHMARK_CASES_DICT
    all_passed = True
    for query, expected_hits in cases.items():
        if query not in drug_features:
            print(f"[{query}] Not in DB!")
            all_passed = False
            continue
            
        res = get_repurposing_candidates(query, top_n=20)
        found = []
        for rank, r in enumerate(res):
            if any(e in r["candidate"] for e in expected_hits):
                found.append((r["candidate"], rank + 1))
                
        if found:
            print(f"[{query}] PASS - Found: {found}")
        else:
            print(f"[{query}] FAIL - No expected hits in Top 20.")
            all_passed = False
            
    return all_passed

def test_hubs():
    print("\n--- Test 3: Hub Assertion ---")
    hubs = ["CLOZAPINE", "OLANZAPINE", "QUETIAPINE", "NEFAZODONE", "FENOFIBRATE", "DICLOFENAC"]
    all_passed = True
    for h in hubs:
        if h not in drug_features:
            print(f"[{h}] Not in DB!")
            all_passed = False
            continue
            
        pct = drug_features[h]["hubness_percentile"]
        if pct > 90.0:
            print(f"[{h}] PASS - Hub Percentile: {pct:.4f}")
        else:
            print(f"[{h}] FAIL - Hub Percentile too low: {pct:.4f}")
            all_passed = False
            
    return all_passed

if __name__ == "__main__":
    t1 = test_chirality()
    t2 = test_benchmark()
    t3 = test_hubs()
    
    if t1 and t2 and t3:
        print("\nALL PHASE 1 GATE CRITERIA PASSED!")
    else:
        print("\nSOME PHASE 1 GATE CRITERIA FAILED!")
