#! python3


__author__ = 'd2jvkpn'
__version__ = '0.2'
__release__ = '2018-09-13'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'


import os

if len(os.sys.argv) == 1 or os.sys.argv[1] in ["-h", "--help"]:
    print ("Extract significantly  expressed genes table, usage:")
    print ("python3 DE_Expression.py <expression.tsv>  <versus.tsv>  <group.tsv>  <outdir>")
    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlicense: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __license__]
    print (_.format (*__))

    os.sys.exit (0)

import pandas as pd

expr, vs, group, outdir = os.sys.argv[1:5]
# expr, vs, group, outdir = "gene_TPM.tsv", "versus.tsv", "group.tsv", "./"
Group = pd.read_csv (group, sep="\t", header=0).to_records()

s2g = {}; g2s = {}

for i in Group:
    _, s, g = i
    if s in s2g:
        s2g[s].append (g)
    else:
        s2g[s] = [g]

    if g in g2s:
        g2s[g].append (s)
    else:
        g2s[g] = [s]

VS = [list(i) for i in pd.read_csv (vs, sep="\t", header=0).values]
Expr = pd.read_csv (expr, sep="\t", index_col=0, header=0)

for i in VS:
    t, u = i
    ss = []
    if t in s2g:
        ss.append (t)
    else:
        ss += g2s[t]

    if u in s2g:
        ss.append (u)
    else:
        ss += g2s[u]

    de = t + "_vs_"+ u + ".DE_signi.tsv"

    tsv = outdir + "/" + t + "_vs_"+ u + ".DE_signi.TPM.tsv"
    gs = pd.read_csv (de, sep="\t", index_col=0).index
    Expr.loc[gs, ss].to_csv (tsv, sep="\t")
