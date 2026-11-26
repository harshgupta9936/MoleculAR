import sqlite3

def fix():
    conn = sqlite3.connect('drugrepo.db')
    c = conn.cursor()
    c.execute("UPDATE drugs SET indication = indication || ', Malaria (Infectious)' WHERE generic_base_name = 'CHLOROQUINE'")
    c.execute("UPDATE drugs SET indication = indication || ', Bacterial Infections (Infectious)' WHERE generic_base_name = 'METRONIDAZOLE'")
    c.execute("UPDATE drugs SET indication = indication || ', Prostate Cancer (Oncology)' WHERE generic_base_name = 'APALUTAMIDE'")
    conn.commit()
    conn.close()
    
    print("Fixed Chloroquine and Metronidazole taxonomy gaps.")

if __name__ == "__main__":
    fix()
