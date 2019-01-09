# pip3 install GEOparse

__author__ = 'd2jvkpn'
__version__ = '0.2'
__release__ = '2019-01-09'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import GEOparse, os
import pandas as pd
from functools import reduce

USAGE = '''Extract GSE(GEO) table, usage:
  python3 GSE_export.py  <gse_id | gse_file>  <output_dir>'''

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']:
    print(USAGE)

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisense: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __license__]

    print (_.format (*__))
    os.sys.exit(2)

####
gse, outdir = os.sys.argv[1], os.sys.argv[2]
if not os.path.isdir(outdir): os.makedirs(outdir)

if gse.endswith(".gz"):
    GSE = GEOparse.get_GEO(filepath = gse)
else:
    GSE = GEOparse.get_GEO(geo=gse, destdir=outdir)

####
tsv = outdir + "/{}.phenotype.tsv".format(GSE.name)
pt = GSE.phenotype_data
pt.to_csv(tsv, sep="\t")
print("Saved", tsv)

ptd, gp = pt["platform_id"].to_dict(), {}

for k in ptd:
    v = ptd[k]
    gp[v] = (gp[v] + [k]) if v in gp else []

####
for k in gp:
    tsv = "{}/{}.{}.infor.tsv".format(outdir, GSE.name, k)
    td = GSE.gpls[k].table

    if td.shape[0] == 0:
        print("no gpls in", GSE.name)
    else:
        td.to_csv(tsv, sep="\t", index=False)
        print("Saved", tsv)

    dfs = []

    for s in gp[k]:
       td = GSE.gsms[s].table

       if td.shape[0] == 0:
           print("not {} expression in {}".format(s, GSE.name))
           continue

       td.iloc[:, 0] = td.iloc[:, 0].astype(str)
       dfs.append(td.iloc[:, [0, 1]])

    if len(dfs) == 0: continue

    DF = reduce(lambda x, y: pd.merge(x, y, 
    how="left", left_on=x.columns[0], right_on=y.columns[0]), dfs)

    tsv = "{}/{}.{}.gsms.tsv".format(outdir, GSE.name, k)
    DF.to_csv(tsv, sep="\t", index=False, na_rep="0")
    print("Saved", tsv)
