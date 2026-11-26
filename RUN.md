# MoleculAR Execution Guide

Quickstart guide for setting up and running the MoleculAR server and dashboard locally.

## 1. Prerequisites & Environment

- **Python:** `3.14.5` (tested on Windows 11 AMD64)
- **Node.js:** `v26.13.0` (with npm)

Create and activate a clean Python virtual environment:

```bash
python -m venv .venv
# On Windows:
.venv\Scripts\activate
# On Linux / macOS:
source .venv/bin/activate

pip install -r requirements.txt
```

## 2. Obtaining the Database

The current SQLite database file `drugrepo.db` is not committed. Earlier builds of `drugrepo.db` exist in the repository history and are NOT the file used for the paper; the SHA-256 identifies the correct file. SIDER is CC BY-NC-SA 4.0 and so the file is not redistributed here.

- **File Name:** `drugrepo.db`
- **Expected SHA-256:** `D88A2333743AF2DBE859B3717764B6A645B9786CAA6F059CC9558E29C7FF0DA7`

The database file is not distributed in this repository. Its availability for reviewers is to be decided; the SHA-256 above identifies the exact file used for the paper.

### Database Build Scripts
`requirements-build.txt` lists the imports of the database-building scripts and has not been install-tested.

### Placement & Verification
Place the database at the repository root:
- **Target Path:** `./drugrepo.db`

Verify the database checksum:
- **Windows (PowerShell):**
  ```powershell
  (Get-FileHash drugrepo.db -Algorithm SHA256).Hash
  ```
- **Linux / macOS:**
  ```bash
  sha256sum drugrepo.db
  ```

## 3. Starting the Server

Start the FastAPI backend service:
```bash
python server.py
```
The server will bind to `http://127.0.0.1:8000`.
Interactive API documentation is available at `http://127.0.0.1:8000/docs`.

## 4. Starting the Frontend Dashboard

In a separate terminal:
```bash
cd frontend
npm install
npm run dev
```
Open your browser to `http://localhost:5173` to access the React dashboard.

### Dashboard Unit Tests
Run the frontend Vitest suite:
```bash
cd frontend
npm test -- --run
```

### Regenerating Benchmark Metrics (`metrics.json`)
The dashboard benchmark panel displays empirical results from `frontend/public/metrics.json`. This file is built directly from the committed duplicate-neutral CSVs (`metrics_C1_P1_dn.csv`, `metrics_C1_P2a_dn_0.4.csv`, `metrics_C1_P3_dn.csv`, `metrics_C1_P4.csv`, and `metrics_C1_P4b.csv`), with the evaluation git tag (`eval-v2`) and target database SHA-256 (`D88A2333...`) displayed in the modal panel.

To regenerate `frontend/public/metrics.json`:
```bash
# Run from repository root:
node frontend/scripts/extract_metrics.js

# Or run from inside frontend/:
cd frontend
node scripts/extract_metrics.js
```

## 5. Running Evaluation & Tests

### API Integration Tests
To test the FastAPI server endpoints, data contracts, and structure rendering:
```bash
pytest tests/test_api.py -v
```

### Evaluation Harness Tests
To execute the evaluation harness test suite:

1. Install the evaluation requirements into your environment:
   ```bash
   pip install -r requirements-eval.txt
   ```

2. The evaluation harness tests in `experiments/tests/` expect the database at `experiments/data/drugrepo_oct6.db`. Copy the database into place:
   ```bash
   # Windows (PowerShell):
   Copy-Item drugrepo.db experiments/data/drugrepo_oct6.db
   # Linux / macOS:
   cp drugrepo.db experiments/data/drugrepo_oct6.db
   ```

3. Run the evaluation pytest suite:
   ```bash
   pytest experiments/tests -v
   ```

## 6. Example API Queries

### Drug Search (Mode A)
Query repurposing candidates for **Sildenafil**:
```bash
curl -X POST http://127.0.0.1:8000/api/rank \
  -H "Content-Type: application/json" \
  -d "{\"drug_name\": \"Sildenafil\", \"top_n\": 5}"
```

### Target Search (Mode B)
Query repurposing candidates for UniProt protein target **P35354** (PTGS2 / COX-2):
```bash
curl -X POST http://127.0.0.1:8000/api/target-search \
  -H "Content-Type: application/json" \
  -d "{\"target_id\": \"P35354\", \"top_n\": 5}"
```
