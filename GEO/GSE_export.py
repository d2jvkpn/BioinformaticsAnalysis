# pip3 install GEOparse

__author__ = 'd2jvkpn'
__version__ = '0.2'
__release__ = '2019-01-09'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import GEOparse, os
import pandas as pd

USAGE = '''Extract GSE(GEO) table, usage:
  python3 GSE_export.py  <gse_id | gse_file>  <output_dir>'''

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']:
    print(USAGE)

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisense: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __license__]

    print (_.format (*__))
    os.sys.exit(2)

def ExprTable (gsms, samples):
    d = gsms[samples[0]].table
    if d.shape[0] == 0: return d

    d = d.iloc[:, [0, 1]]
    idxn = d.columns[0]

    for i in samples[1:]:
        t = GSE.gsms[i].table.iloc[:, [0, 1]]
        t.iloc[:, [0]] = t.iloc[:, 0].astype(str)
        h = list(t.columns); h = [idxn, i]; t.columns = h
        d = pd.merge(d, t, how="left", left_on=idxn, right_on=idxn)

    return d

gse, outdir = os.sys.argv[1], os.sys.argv[2]
if not os.path.isdir(outdir): os.makedirs(outdir)

if gse.endswith(".gz"):
    GSE = GEOparse.get_GEO(filepath = gse)
else:
    GSE = GEOparse.get_GEO(geo=gse, destdir=outdir)

tsv = outdir + "/{}.phenotype.tsv".format(GSE.name)
pt = GSE.phenotype_data
pt.to_csv(tsv, sep="\t")
print("Saved", tsv)

ED = ExprTable(GSE.gsms, list(GSE.gsms))

if ED.shape[0] == 0:
    os.sys.exit("Not GSMS available in " + gse)

ptd, gp = pt["platform_id"].to_dict(), {}

for k in ptd:
    v = ptd[k]
    gp[v] = (gp[v] + [k]) if v in gp else []

for k in gp:
    tsv = "{}/{}.{}.infor.tsv".format(outdir, GSE.name, k)
    GSE.gpls[k].table.to_csv(tsv, sep="\t", index=False)
    print("Saved", tsv)

    tsv = "{}/{}.{}.gsms.tsv".format(outdir, GSE.name, k)
    ED.loc[:, [ED.columns[0]] + gp[k]].to_csv(tsv, sep="\t", index=False)
    print("Saved", tsv)
