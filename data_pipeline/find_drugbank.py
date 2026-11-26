import os
import sqlite3

CHEMBL_DB_PATH = os.getenv("CHEMBL_DB_PATH", os.path.join(os.path.dirname(__file__), "..", "chembl_37", "chembl_37_sqlite", "chembl_37.db"))
conn = sqlite3.connect(CHEMBL_DB_PATH)
cursor = conn.cursor()

cursor.execute("SELECT syn_type, COUNT(1) FROM molecule_synonyms GROUP BY syn_type")
for row in cursor.fetchall():
    print(row)

# also check cross-references table if it exists
cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
tables = [r[0] for r in cursor.fetchall()]
if 'chembl_id_lookup' in tables:
    print("chembl_id_lookup exists")
if 'compound_records' in tables:
    print("compound_records exists")
    cursor.execute("SELECT src_id, COUNT(1) FROM compound_records GROUP BY src_id")
    print(cursor.fetchall())
conn.close()
