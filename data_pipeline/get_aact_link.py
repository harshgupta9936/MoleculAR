import urllib.request, re
url = 'https://aact.ctti-clinicaltrials.org/download'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as r:
        content = r.read().decode('utf-8')
        links = re.findall(r'href=[\'\"]?([^\'\" >]+_pipe-delimited-export\.zip)', content)
        if links:
            link = links[0]
            if not link.startswith('http'):
                link = 'https://aact.ctti-clinicaltrials.org' + link
            print('Found link:', link)
        else:
            print('No link found')
except Exception as e:
    print('Failed:', e)
