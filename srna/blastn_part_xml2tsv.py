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

            results.append([r.query, a.hit_def, ident, h.align_length,
            h.match.count(" "), h.gaps, h.query_start, h.query_end, 
            h.sbjct_start, h.sbjct_end, h.expect, h.score, h.match])

# a.hit_id, h.query

f.close()

d1 = pd.DataFrame(results)
d1.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
"qstart", "qend", "sstart", "send", "evalue", "bitscore", "midline"]

d1["qseqid"] = d1["qseqid"].apply(lambda x: x.split()[0])
d1 = d1[d1["gapopen"] == 0]

def parseMidline(s):
    m = ""
    for i in s.split():
       if len(i) > len(m): m = i

    return (s.index(m), len(m))

d1["perfect"] = d1["midline"].apply(lambda x: parseMidline(x))

def newAlign(m):
    if m["mismatch"] == 0:
        return (m["length"], m["qstart"], m["qend"], m["sstart"], m["send"])

    length = m["perfect"][1]
    qstart = m["qstart"] + m["perfect"][0]
    qend = qstart + length - 1

    if m["sstart"] < m["send"]:
         sstart = m["sstart"] + m["perfect"][0]
         send = sstart + length - 1
    else:
        sstart = m["sstart"] - m["perfect"][0]
        send = sstart - length + 1

    return (length, qstart, qend, sstart, send)


d1["new_align"] = d1.apply(lambda x: newAlign(x), axis=1)

d1["length"] = d1["new_align"].apply(lambda x: x[0])
d1["qstart"] = d1["new_align"].apply(lambda x: x[1])
d1["qend"] = d1["new_align"].apply(lambda x: x[2])
d1["sstart"] = d1["new_align"].apply(lambda x: x[3])
d1["send"] = d1["new_align"].apply(lambda x: x[4])

d1["pident"], d1["mismatch"] = 100, 0

os.makedirs(os.path.dirname(os.path.abspath(tsv)), exist_ok=True)
d1.drop(["midline", "perfect", "new_align"], axis=1).to_csv(tsv, sep="\t", index=False)
