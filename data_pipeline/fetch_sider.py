import pandas as pd
import urllib.request
import os
import json

def fetch_sider():
    files = {
        "meddra_all_se.tsv.gz": "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz",
        "drug_names.tsv": "http://sideeffects.embl.de/media/download/drug_names.tsv"
    }
    
    for filename, url in files.items():
        if not os.path.exists(filename):
            print(f"Downloading {filename}...")
            urllib.request.urlretrieve(url, filename)

    print("Loading SIDER data... (this might take a few seconds)")
    sider = pd.read_csv("meddra_all_se.tsv.gz", sep="\t", header=None,
                        names=["stitch_id_flat", "stitch_id_stereo", "umls_cid",
                               "meddra_type", "umls_cid_meddra", "side_effect"])
    
    drug_names = pd.read_csv("drug_names.tsv", sep="\t", header=None,
                             names=["stitch_id", "drug_name"])

    print("Building side-effect sets...")
    drug_side_effects = sider.groupby("stitch_id_flat")["side_effect"].apply(set).to_dict()

    name_to_stitch = {}
    for _, row in drug_names.iterrows():
        name_to_stitch[str(row["drug_name"]).upper()] = row["stitch_id"]

    output_se = {}
    for name, stitch_id in name_to_stitch.items():
        if stitch_id in drug_side_effects:
            output_se[name] = list(drug_side_effects[stitch_id])
            
    with open("drug_side_effects.json", "w") as f:
        json.dump(output_se, f, indent=2)

    print(f"Saved side-effect profiles for {len(output_se)} named drugs to drug_side_effects.json")

if __name__ == "__main__":
    fetch_sider()
