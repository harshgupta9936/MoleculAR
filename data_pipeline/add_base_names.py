import sqlite3

SUFFIXES = [
    " HYDROCHLORIDE", " SODIUM", " POTASSIUM", " CALCIUM", " MAGNESIUM", " ZINC",
    " SULFATE", " MESYLATE", " BESYLATE", " TOSYLATE", " MALEATE", " FUMARATE",
    " TARTRATE", " CITRATE", " ACETATE", " PHOSPHATE", " DIPHOSPHATE", " NITRATE",
    " BROMIDE", " HYDROBROMIDE", " CHLORIDE", " IODIDE", " FLUORIDE", " LACTATE", " SUCCINATE",
    " PAMOATE", " VALERATE", " STEARATE", " PALMITATE", " OLEATE", " SALICYLATE",
    " BENZOATE", " NAPSYLATE", " GLUCONATE", " GLUCEPTATE", " LACTOBIONATE",
    " ANHYDROUS", " MONOHYDRATE", " DIHYDRATE", " TRIHYDRATE", " TETRAHYDRATE",
    " PENTA", " HEXA", " HEPTA", " OCTA", " NAI", " FREE ACID", " FREE BASE",
    " CATION", " ANION", ", (E)-", ", (Z)-"
]

def clean_drug_name(name):
    cleaned = str(name).upper()
    changed = True
    while changed:
        changed = False
        for suffix in SUFFIXES:
            if cleaned.endswith(suffix):
                cleaned = cleaned[:-len(suffix)].strip()
                changed = True
    return cleaned

conn = sqlite3.connect("drugrepo.db")
cursor = conn.cursor()

# Check if column exists
cursor.execute("PRAGMA table_info(drugs)")
columns = [col[1] for col in cursor.fetchall()]

if "generic_base_name" not in columns:
    cursor.execute("ALTER TABLE drugs ADD COLUMN generic_base_name VARCHAR(255)")

cursor.execute("SELECT chembl_id, drug_name FROM drugs")
rows = cursor.fetchall()

print(f"Populating generic_base_name for {len(rows)} drugs...")
for cid, name in rows:
    base = clean_drug_name(name)
    cursor.execute("UPDATE drugs SET generic_base_name = ? WHERE chembl_id = ?", (base, cid))

conn.commit()
conn.close()
print("Done!")
