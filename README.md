# MoleculAR: AI-Powered Drug Repurposing Screener

A computational framework for screening repurposing candidates among FDA-approved small-molecule drugs using multi-signal similarity scoring and machine learning.

## System Architecture & Tech Stack

MoleculAR is architected as a decoupled client-server web application:
- **Backend API:** FastAPI (`server.py`) running on Python 3.11+, serving RESTful endpoints for candidate screening, target querying, and validation metrics on `http://localhost:8000`.
- **Frontend Dashboard:** React single-page application built with Vite (`frontend/`), featuring interactive candidate search, structural visualization, and feature explainability panels on `http://localhost:5173`.
- **Computational Core:** RDKit (ECFP4 Morgan fingerprinting, SMILES parsing), Scikit-Learn, and XGBoost.
- **Data Layer:** SQLite database (`drugrepo.db`).

*(Note: Earlier project notes proposed a Streamlit interface via `app.py`, but the production implementation is entirely built with FastAPI and React).*

## Pipeline & Data Sources

Candidates are ranked using a multi-signal scoring engine combining three pharmacological signals:
1. **Chemical Structure Similarity:** Morgan Fingerprints (ECFP4, 2048-bit, radius 2) compared via Tanimoto coefficient (weight = 0.4).
2. **Biological Target Overlap:** Jaccard similarity over shared human protein targets mapped to UniProt IDs (weight = 0.4).
3. **Phenotypic Side-Effect Similarity:** Jaccard similarity over shared adverse event terms mapped to MedDRA vocabulary (weight = 0.2).
4. **Hub Bias Normalization:** Percentile penalty suppressing promiscuous structural hubs based on average Morgan fingerprint Tanimoto centrality.

### Data Sources & Coverage
- **ChEMBL:** 2,267 FDA-approved small-molecule drugs with validated SMILES structures and biological targets.
- **SIDER:** Adverse drug reactions mapped to MedDRA vocabulary.
*(Exploratory proposals for DrugBank, CTD, and DisGeNET were not incorporated into the operational scoring pipeline).*

### Database Quality & Known Limitations
- **Trivial Compounds:** Exactly 65 trivial compounds (25 radioisotopes, 22 inorganic salts, 9 medical gases, 9 tiny organics) are tagged; 64 are excluded from structural similarity scoring (leaving 2,202 scored drugs). Sodium Chloride (`CHEMBL1200574`) is noted as mis-tagged as Tiny Organic.
- **Zero-Target Coverage:** 693 drugs (30.57%) have zero documented biological targets in `drugrepo.db`.
- **Indication & ATC Gaps:** Therapeutic indication and ATC classification are recorded as "Other" for 91.4% of compounds.
- **Benchmark Data Gaps:** Levofloxacin (`CHEMBL33`) and Pinacidil (`CHEMBL1159`) have 0 targets and 0 side effects in the database.

## Measured Benchmark Results

All evaluation metrics are verified against reproducible test scripts in `/verification/`:
- **Retrospective Benchmark (8 Canonical Cases):**
  - **Hand-Weighted Heuristic:** **6/8 (75.0%) hit rate at Top-10 and Top-20**.
    - Hits: Sildenafil → Vardenafil (#1), Minoxidil → Diazoxide (#1), Metformin → Phenformin (#9), Thalidomide → Lenalidomide (#1), Duloxetine → Fluoxetine (#3), Gabapentin → Pregabalin (#2).
    - Misses: Bupropion → Varenicline (Rank 66), Raloxifene → Fulvestrant (Rank >100; Dienestrol Rank >100).
  - **Supervised ML (XGBoost):** **5/8 (62.5%) hit rate at Top-10** (misses Metformin #36, Bupropion #230, Raloxifene #216).
- **Cross-Validation Performance (StratifiedGroupKFold, 5 folds grouped by `drug_name` across 156 unique drugs):**
  - **XGBoost:** ROC-AUC **0.8564 ± 0.0622**
  - **Logistic Regression:** ROC-AUC **0.8615 ± 0.0291**
  - **Unfitted Heuristic:** ROC-AUC **0.8061 ± 0.0343**
  *(Note: The previously claimed 0.9392 XGBoost AUC was inflated due to row-level data leakage across identical drugs; true grouped generalization is ~0.856–0.862).*

## Setup & Execution

### 1. Environment Setup
```bash
# Activate virtual environment
.\drugrepo-env\Scripts\activate

# Install backend dependencies
pip install -r requirements.txt
```

### 2. Start Backend API
```bash
python server.py
# API server starts at http://127.0.0.1:8000
# OpenAPI documentation at http://127.0.0.1:8000/docs
```

### 3. Start Frontend Dashboard
```bash
cd frontend
npm install
npm run dev
# Vite dev server runs at http://localhost:5173
```

### 4. Run Verification Suite
```bash
python data_pipeline/run_validation.py
python train_and_eval_xgboost_v2.py
python verification/step3_verification.py
```

## Reproducibility

* **DB-A Artefact:** The primary database artefact DB-A is not stored in Git. Its SHA-256 checksum is `D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7` (recorded in `experiments/results/manifest.json`). The database artefact is available on request or to be deposited in a research data repository.
* **Data Pipeline:** The scripts in `data_pipeline/` document how the dataset was built and serve as provenance documentation, not a guarantee of a bit-identical rebuild.
* **Evaluation Suite:** Retrospective benchmark evaluation and protocols live in `experiments/` (see `experiments/REPRODUCE.md`).
* **Scoring Engine Defaults & Hub Penalty:** The similarity weights ($w_s = 0.4, w_t = 0.4, w_e = 0.2$) in `scoring_engine.py` and the evaluated production hub penalty $1 + (\text{hubness\_percentile} / 100)^2$ in `pipeline.py` (evaluated as Method M7 in `eval_engine.py`) are the production defaults; retrospective evaluation across protocols shows no consistent benefit of the hub penalty (see `experiments/RESULTS.md`).

## Disclaimer
MoleculAR is an in-silico computational screening tool for hypothesis generation. It does not provide medical advice or guarantee clinical efficacy. All predicted candidates require rigorous preclinical in vitro assays and clinical trial validation.
## Data and licences

The scored database (`drugrepo.db`) is **not distributed** in this repository. Its expected SHA-256 checksum is:

`D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7`

Dataset provenance, versions, and access metadata are documented in [`experiments/DATA_SOURCES.md`](experiments/DATA_SOURCES.md). ChEMBL 37 is licensed under CC BY-SA 3.0, DrugCentral under CC BY-SA 4.0, and SIDER 4.1 under CC BY-NC-SA 4.0, so the database built from them is not redistributed here.

The PostgreSQL credentials referenced in the evaluation and data-extraction scripts (`host="unmtid-dbs.net"`, `port=5433`, `user="drugman"`) are DrugCentral's own public read-only login.

| File Group | Paths | Licence & Attribution |
| :--- | :--- | :--- |
| **Code** | Python backend, evaluation harness, data-pipeline scripts, tests, and React frontend (`*.py`, `*.js`, `*.jsx`, `*.css`, `*.html`, `*.sh`, `*.sql`) | **MIT** |
| **DrugCentral-derived files** | Raw extracts (`experiments/data/drugcentral_raw_approved_indications.csv`, `drugcentral_raw_atc.csv`), `gold_source.txt`, gold-pair and mapping CSVs (`gold_pairs*.csv`, `excluded_indications.csv`, `rejected_rule_b.csv`, `rule_b_matches.csv`, `rule_c_matches.csv`, `many_to_one.csv`, `sampled_20_with_names.csv`), and `experiments/data/unmapped.csv` | **CC BY-SA 4.0** — attribution to DrugCentral (Avram et al., *Nucleic Acids Res.* 2023) |
| **SIDER-derived files** | `experiments/results/se_leakage_terms.csv`, `experiments/results/se_top50_lift.csv` | **CC BY-NC-SA 4.0** — attribution to SIDER 4.1 |
| **Result CSVs (numeric metrics only)** | Aggregate evaluation tables (`experiments/results/metrics_*.csv`, `coverage_*.csv`, `paired_differences*.csv`, `hub_penalty_ablation.csv`, `weight_*.csv`, `tie_sensitivity*.csv`, `shuffle_control_results.csv`, `se_leakage_ablation.csv`, `se_leakage_approved_ablations.csv`, `p4_*query_census.csv`, `known_cases_*.csv`, `failure_analysis_*.csv`, `duplicate_pairs.csv`) and `frontend/public/metrics.json` | **MIT** |
| **Figures and screenshots** | `experiments/figures/*.png`, `screenshots/*.png` | **MIT** |

### Tags

The tags `eval-v1`, `eval-v2`, and `demo-v1` were recreated on snapshot commits of a private development history. The corresponding source commit hashes (`a732c9c`, `2aec56c`, `c9ff106`, and `a4708a2`) are listed in the commit messages, and result integrity is verified by the SHA-256 values in [`experiments/results/manifest.json`](experiments/results/manifest.json).
