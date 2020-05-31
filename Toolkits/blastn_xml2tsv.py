# https://github.com/d2jvkpn/BioinformaticsAnalysis

import os
from Bio.Blast import NCBIXML
import pandas as pd

xml, tsv = os.sys.argv[1:3]
# xml, prefix = "Aa.random_nt.xml", "Aa.blast.tsv"

results = []

f = open(xml, "r")

for r in NCBIXML.parse(f):
    for a in r.alignments:
        for h in a.hsps:
            ident = "%.2f%%" % (h.identities/min(len(h.query), len(h.sbjct))*100)

            results.append([r.query, a.hit_id, ident, h.align_length,
            h.match.count(" "), h.gaps, h.query_start, h.query_end, 
            h.sbjct_start, h.sbjct_end, h.expect, h.score, a.hit_def, h.query, h.match])

# a.hit_id, h.query

f.close()

d1 = pd.DataFrame(results)
d1.columns = ["query", "hit_id", "ident", "align_length", "missmatch", "gaps",
"query_start", "query_end", "sbjct_start", "sbjct_end", "evalue", "bitscore", 
"hit_def", "qseq", "midline"]

os.makedirs(os.path.dirname(os.path.abspath(tsv)), exist_ok=True)
d1.to_csv(tsv, sep="\t", index=False)
