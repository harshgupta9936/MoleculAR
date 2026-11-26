import ftplib
import os
import gzip

def download_file(ftp, remote_path, local_path):
    print(f"Downloading {remote_path}...")
    with open(local_path, 'wb') as f:
        ftp.retrbinary(f'RETR {remote_path}', f.write)

try:
    ftp = ftplib.FTP('ftp.ebi.ac.uk')
    ftp.login()
    
    os.makedirs('unichem_data', exist_ok=True)
    
    # Download ChEMBL to DrugBank (src1 to src2)
    download_file(ftp, '/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src2.txt.gz', 'unichem_data/src1src2.txt.gz')
    
    # Download ChEMBL to PubChem (src1 to src22)
    # Wait, let's verify pubchem is src 22
    download_file(ftp, '/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz', 'unichem_data/src1src22.txt.gz')
    
    ftp.quit()
    print("Download complete.")
except Exception as e:
    print("FTP error:", e)
