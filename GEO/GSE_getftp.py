__author__ = 'd2jvkpn'
__version__ = '0.2'
__release__ = '2019-01-08'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import GEOparse, os
import pandas as pd

USAGE = '''Scrap GSE(GEO) family.soft.gz ftp link, usage:
  python3  GSE_getftp.py  <gse_id...>'''

if len(os.sys.argv) ==1 or os.sys.argv[1] in ['-h', '--help']:
    print(USAGE)

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisense: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __license__]

    print (_.format (*__))

    os.sys.exit(2)

import os, requests, sys
from bs4 import BeautifulSoup

def GetLink (gse):
    links = []
    prefix = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
    query = requests.get(prefix + gse)

    if query.status_code != 200:
        sys.stderr.write("failed to get GSE(%s) page" % gse)
        return links

    bs = BeautifulSoup(query.text, 'html.parser')

    lk = bs.find("a", string="SOFT formatted family file(s)").get("href")
    links.append(lk + "/{}_family.soft.gz".format(gse))

    trs = bs.find("td", string="Supplementary file").parent.parent.find_all("tr")

    for tr in trs[1:]:
        tds = tr.find_all("td")
        if len(tds) < 4: continue
        lk = tds[2].find("a").get("href")
        if lk.startswith("ftp:"): links.append(lk)

    return links

for gse in os.sys.argv[1:]:
    links = GetLink(gse)
    for lk in links: print("{}  {}".format(gse, lk))
