import os, time

__author__ = 'd2jvkpn'
__version__ = '0.3'
__release__ = '2019-01-26'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

USAGE = '''Scraping circRNA-miRNA target from circinteractome and save the result to tsv
Usage:
  python3  circinteractome_miRNA.py  <circRNA_list_file>  <outdir>
  Note: pandas, selenium, firefox and geckodriver are required 
'''

if len(os.sys.argv) == 1 or os.sys.argv[1] in ['-h', '--help']:
    print(USAGE)

    _ = 'author:  {}\nversion: {}\nrelease: {}\nproject: {}\nlicense: {}'
    __ = [__author__,  __version__, __release__, __project__, __license__]
    print (_.format (*__))
    os.sys.exit(2)

####
import pandas as pd
from selenium import webdriver

URL = 'https://circinteractome.nia.nih.gov/miRNA_Target_Sites/mirna_target_sites.html'

circlist, outdir = os.sys.argv[1], os.sys.argv[2]

with open(circlist) as f: circs = f.read().splitlines()

if not os.path.isdir(outdir): os.makedirs(outdir)

####
options = webdriver.FirefoxOptions()
options.set_headless()    # options.add_argument('-headless')
options.add_argument('--disable-gpu')
firefox = webdriver.Firefox(firefox_options=options)

def Arich(firefox, URL, c, tsv):
    print("Querying {}, {}".format(c, time.strftime("%Y-%m-%d %H:%M:%S %z")))
    firefox.get(URL)
    firefox.find_element_by_name("gcircrna").send_keys(c)
    firefox.find_element_by_name("submit").click()

    tbls = pd.read_html(firefox.page_source)
    # <table> block is nesting (not standard)
    nr = [tb.shape[0] for tb in tbls]
    tbl = tbls[nr.index(max(nr))]

    tbl = tbl.iloc[:, 0:11]
    tbl = tbl.loc[tbl.iloc[:, 0].dropna().index, :]
    header = list(tbl.loc[tbl.iloc[:, 0] == "CircRNAMirbase ID", :].iloc[0, :])
    index = [i.startswith("hsa_circ_") for i in tbl.iloc[:, 0]]
    tbl = tbl.loc[index, :]
    tbl.columns = header
    tbl.iloc[:, 0] = [i.replace(u'\xa0', u' ') for i in tbl.iloc[:, 0]]

    tbl.to_csv(tsv, sep="\t", index=False)
    print("    saved {}, {} records".format(tsv, tbl.shape[0]))

for c in circs:
    tsv = "{}/circinteractome_miRNA__{}.tsv".format(outdir, c)
    if os.path.isfile(tsv) or not c.startswith("hsa_circ_"):
        continue

    try:
        Arich(firefox, URL, c, tsv)
    except:
        print("    failed to achive", c)

    # time.sleep(5)

firefox.close()
