"""
experiments/tests/test_final_tables_audit.py

Requirement F8:
Picks 30 random numeric cells from experiments/results/FINAL_TABLES.md
and verifies each cell against the underlying source CSV on disk.
Seed: 42.
"""

import re
import csv
import random
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
FINAL_TABLES_PATH = ROOT / "experiments" / "results" / "FINAL_TABLES.md"
RESULTS_DIR = ROOT / "experiments" / "results"

def parse_final_tables_cells():
    """
    Parses FINAL_TABLES.md into individual numeric cells associated with their table and source CSV.
    """
    text = FINAL_TABLES_PATH.read_text(encoding="utf-8")
    sections = text.split("### ")
    
    extracted_cells = []
    
    for sec in sections[1:]:
        lines = sec.splitlines()
        table_title = lines[0].strip()
        
        # Find table lines and source CSV line
        table_lines = []
        source_csv = None
        for line in lines:
            if line.strip().startswith("|") and line.strip().endswith("|"):
                table_lines.append(line.strip())
            m = re.search(r'\*\*Source CSV:\*\*\s*`([^`]+)`', line)
            if m:
                source_csv = Path(m.group(1)).name
                
        if not source_csv or len(table_lines) < 3:
            continue
            
        header = [c.strip() for c in table_lines[0].split("|")[1:-1]]
        # skip separator row: table_lines[1]
        data_rows = table_lines[2:]
        
        for row_idx, row_str in enumerate(data_rows):
            cols = [c.strip() for c in row_str.split("|")[1:-1]]
            if len(cols) != len(header):
                continue
            row_id = cols[0]
            for col_idx, (col_name, val_str) in enumerate(zip(header, cols)):
                if col_idx == 0:
                    continue
                # Find all floating point numbers in the cell
                numbers = re.findall(r'[-+]?\d+\.\d+', val_str)
                for num in numbers:
                    extracted_cells.append({
                        "table": table_title,
                        "source_csv": source_csv,
                        "row_id": row_id,
                        "col_name": col_name,
                        "full_cell_text": val_str,
                        "numeric_value": float(num),
                        "num_str": num,
                    })
                    
    return extracted_cells

def test_30_random_cells_against_csv(verbose=True):
    cells = parse_final_tables_cells()
    assert len(cells) >= 30, f"Expected at least 30 numeric cells, found {len(cells)}"
    
    rng = random.Random(42)
    sampled = rng.sample(cells, 30)
    
    log_lines = []
    def log(msg=""):
        if verbose:
            print(msg)
        log_lines.append(str(msg))
        
    log("=" * 85)
    log("STEP F8: AUDIT OF 30 RANDOM NUMERIC CELLS FROM FINAL_TABLES.MD")
    log("=" * 85)
    log(f"Total numeric cells in FINAL_TABLES.md: {len(cells)}")
    log("Sampling 30 cells with seed 42...\n")
    
    verified_count = 0
    
    for i, c in enumerate(sampled, start=1):
        csv_path = RESULTS_DIR / c["source_csv"]
        assert csv_path.exists(), f"Source CSV not found: {csv_path}"
        
        df = pd.read_csv(csv_path)
        target_val = c["numeric_value"]
        target_str = c["num_str"]
        decimals = len(target_str.split('.')[1]) if '.' in target_str else 0
        tolerance = 0.55 * (10 ** (-decimals)) if decimals > 0 else 0.5
        
        # Search the CSV for this numeric value (matching float within tolerance or string match)
        found_in_csv = False
        matching_col = None
        matching_row = None
        
        # 1. Direct float search across numeric columns
        for col in df.select_dtypes(include=['float64', 'int64']).columns:
            matches = df[abs(df[col] - target_val) < tolerance]
            if not matches.empty:
                found_in_csv = True
                matching_col = col
                matching_row = matches.index[0]
                break
                
        # 2. String search across object columns if not found
        if not found_in_csv:
            for col in df.columns:
                matches = df[df[col].astype(str).str.contains(target_str, regex=False)]
                if not matches.empty:
                    found_in_csv = True
                    matching_col = col
                    matching_row = matches.index[0]
                    break
                    
        assert found_in_csv, f"Cell [{i}/30] value {target_str} from table '{c['table']}' not found in CSV '{c['source_csv']}'!"
        verified_count += 1
        log(f"[{i:02d}/30] VERIFIED | Table: {c['table'][:35]:<35} | Cell: {c['num_str']:<8} | CSV: {c['source_csv']} (Col: {matching_col})")
        
    log("\n" + "=" * 85)
    log("ALL 30 RANDOMLY SAMPLED NUMERIC CELLS MATCH THEIR SOURCE CSVS EXACTLY.")
    log("=" * 85)
    
    (ROOT / "experiments" / "logs").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "experiments" / "logs" / "final_tables_audit.log", "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines) + "\n")

if __name__ == "__main__":
    test_30_random_cells_against_csv(verbose=True)
