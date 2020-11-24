import os, sys, time

__author__ = 'd2jvkpn'
__version__ = '0.5'
__release__ = '2019-02-17'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

USAGE = '''Scraping circRNA-miRNA target from circinteractome and save the result to tsv
Usage:
  $ python3  biocrawler_circinteractome_miRNA.py  <circRNA_list_file>  <outdir>
  Note: pandas, selenium, geckodriver and firefox are required 
'''

if len(os.sys.argv) < 3 or os.sys.argv[1] in ['-h', '--help']:
    print(USAGE)

    _ = 'author:  {}\nversion: {}\nrelease: {}\nproject: {}\nlicense: {}'
    __ = [__author__,  __version__, __release__, __project__, __license__]
    print (_.format (*__))
    os.sys.exit(2)

####
import pandas as pd
from selenium import webdriver
from bs4 import BeautifulSoup
from signal import signal, SIGINT, SIG_DFL
signal(SIGINT, SIG_DFL)

URL = 'https://circinteractome.nia.nih.gov/miRNA_Target_Sites/mirna_target_sites.html'

circlist, outdir = os.sys.argv[1], os.sys.argv[2]

with open(circlist) as f: circs = f.read().splitlines()

if not os.path.isdir(outdir): os.makedirs(outdir)

####
options = webdriver.FirefoxOptions()
options.set_headless()    # options.add_argument('-headless')
options.add_argument('--disable-gpu')
firefox = webdriver.Firefox(firefox_options=options)

def Query(firefox, URL, c):
    firefox.get(URL)
    firefox.find_element_by_name("gcircrna").send_keys(c)
    firefox.find_element_by_name("submit").click()

    hs = BeautifulSoup(firefox.page_source, 'html.parser')
    tbl = pd.read_html(str(hs.find('table', width='1200')))[0]

    n1 = list(tbl.iloc[:, 0] == "CircRNAMirbase ID").index(True)
    tbl.columns = tbl.iloc[n1, :]
    tbl = tbl.iloc[(n1+1):, :]
    tbl.iloc[:, 0] = [i.replace(u'\xa0', u' ') for i in tbl.iloc[:, 0]]
    return tbl


for c in circs:
    tsv = "{}/circinteractome_miRNA__{}.tsv".format(outdir, c)
    if os.path.isfile(tsv) or not c.startswith("hsa_circ_"):
        continue

    msg = "Querying {}, {}".format(c, time.strftime("%Y-%m-%d %H:%M:%S %z"))
    print(msg, file=sys.stdout, flush=True)

    try:
        tbl = Query(firefox, URL, c)
        tbl.to_csv(tsv, sep="\t", index=False)
        msg = "saved results of {}, {} records".format(c, tbl.shape[0])
        print("    " + msg, file=sys.stdout, flush=True)
    except:
        msg = "!!! failed to achive {}".format(c)
        print("    " + msg, file=sys.stderr, flush=True)

firefox.close()
