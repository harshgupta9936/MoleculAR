import sqlite3
import json
from pipeline import get_repurposing_candidates, init_db_cache, drug_features

def print_hubs():
    print("--- NEW HUB LIST (Top 5 > 90th Percentile) ---")
    conn = sqlite3.connect("drugrepo.db")
    cursor = conn.cursor()
    cursor.execute("SELECT generic_base_name, hubness_percentile FROM drugs ORDER BY hubness_percentile DESC LIMIT 10")
    for row in cursor.fetchall():
        if row[1] > 90.0:
            print(f"{row[0]}: {row[1]:.2f}")
    conn.close()

def rank_and_find(query, expected_hits):
    print(f"\n--- {query} TOP 20 & EXPECTED HITS ---")
    res = get_repurposing_candidates(query, top_n=3000) # Get all to find rank
    
    print("Top 20:")
    for rank, r in enumerate(res[:20]):
        print(f"  {rank+1}. {r['candidate']} (Score: {r['score']:.4f}) | " 
              f"Str: {r['struct_sim']:.4f}, "
              f"Tar: {r['target_sim']:.4f}, "
              f"SE: {r['side_effect_sim']:.4f}")
              
    print("\nExpected Hits:")
    for expected in expected_hits:
        found = False
        for rank, r in enumerate(res):
            if expected in r['candidate']:
                print(f"  {expected} -> Rank {rank+1} | Score: {r['score']:.4f} "
                      f"(Str: {r['struct_sim']:.4f}, Tar: {r['target_sim']:.4f}, SE: {r['side_effect_sim']:.4f})")
                found = True
                break
        if not found:
            print(f"  {expected} -> NOT FOUND IN RESULTS")

if __name__ == "__main__":
    init_db_cache()
    print_hubs()
    rank_and_find("BUPROPION", ["VARENICLINE"])
    rank_and_find("RALOXIFENE", ["FULVESTRANT", "DIENESTROL"])
