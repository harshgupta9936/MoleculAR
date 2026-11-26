import urllib.request
import datetime

for days_back in range(5):
    d = datetime.date.today() - datetime.timedelta(days=days_back)
    date_str = d.strftime("%Y%m%d")
    url = f"https://aact.ctti-clinicaltrials.org/static/exported_files/daily/{date_str}_pipe-delimited-export.zip"
    print("Trying", url)
    try:
        req = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(req) as r:
            if r.status == 200:
                print("Found:", url)
                break
    except Exception as e:
        print("Error:", e)
