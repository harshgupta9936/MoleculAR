import urllib.request
import zipfile
import io
import os

url = "https://aact.ctti-clinicaltrials.org/static/exported_files/daily/20261003_pipe-delimited-export.zip"
print(f"Downloading AACT bulk dataset from {url}...")

# We will stream the download and write to a local zip file
zip_path = "aact_export.zip"
if not os.path.exists(zip_path):
    urllib.request.urlretrieve(url, zip_path)

print("Download complete. Extracting specific files...")
with zipfile.ZipFile(zip_path, 'r') as z:
    for filename in ['interventions.txt', 'conditions.txt']:
        if filename in z.namelist():
            print(f"Extracting {filename}...")
            z.extract(filename, path='.')
            
print("Extraction complete.")
