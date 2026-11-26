import urllib.request, re
url = 'https://aact.ctti-clinicaltrials.org/'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as r:
        content = r.read().decode('utf-8')
        links = re.findall(r'href=[\"\']([^\"\']+)[\"\']', content)
        print('Links found:', [l for l in links if 'download' in l or 'pipe' in l or 'zip' in l])
except Exception as e:
    print('Failed:', e)
