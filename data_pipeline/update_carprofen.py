import sqlite3

conn = sqlite3.connect('drugrepo.db')
cursor = conn.cursor()

# Find Carprofen
cursor.execute("SELECT chembl_id, drug_name, market_status FROM drugs WHERE generic_base_name='CARPROFEN'")
print('Before update:', cursor.fetchall())

cursor.execute("UPDATE drugs SET market_status='Manual Review - Veterinary (Non-US)' WHERE generic_base_name='CARPROFEN'")
conn.commit()

cursor.execute("SELECT chembl_id, drug_name, market_status FROM drugs WHERE generic_base_name='CARPROFEN'")
print('After update:', cursor.fetchall())
conn.close()
