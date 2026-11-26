import sqlite3

def fix():
    conn = sqlite3.connect('drugrepo.db')
    c = conn.cursor()
    c.execute("UPDATE drugs SET indication = 'Prostate Cancer (Oncology)' WHERE generic_base_name = 'APALUTAMIDE'")
    c.execute("UPDATE drugs SET indication = 'Hepatitis C (Infectious), Cirrhosis (Gastroenterology)' WHERE generic_base_name = 'DACLATASVIR'")
    c.execute("UPDATE drugs SET indication = 'Gastroenteropancreatic neuroendocrine tumors (Oncology, Gastroenterology)' WHERE generic_base_name = 'LANREOTIDE'")
    conn.commit()
    
    c.execute("SELECT chembl_id, generic_base_name, indication FROM drugs WHERE indication != 'Other' AND indication IS NOT NULL AND indication != '' AND trivial_category IS NULL")
    drugs_raw = c.fetchall()
    conn.close()
    
    print("Fixed Apalutamide, Daclatasvir, and Lanreotide.")

if __name__ == "__main__":
    fix()
