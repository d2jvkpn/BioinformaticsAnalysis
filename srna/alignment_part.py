import os
import pandas as pd
from Bio import SeqIO

## glob
## fs = glob.glob("blastn/*.raw.tsv")
minLen = int(os.sys.argv[1])
fs = os.sys.argv[2:]

for f in fs:
    s = os.path.basename(f).replace(".raw.tsv.gz", "").replace("blastn_", "")
    print("Processing alignment of", s)
    out = f.replace(".raw.tsv.gz", ".align.tsv")

    bd = pd.read_csv(f, sep="\t")
    bd = bd[bd["mismatch"] == 0]
    bd = bd[bd["pident"] == 100] # no mismatch and gapopen in aligned fragment 
    bd = bd[bd["length"] >= minLen]
    h = list(bd.columns); h[0] = "id"; bd.columns = h
    idx = bd["id"].drop_duplicates().index
    bd = bd.loc[idx, :]

    sd = pd.read_csv("clean/%s.fasta.tsv.gz" % s, sep="\t").drop(["length"], axis=1)
    h = list(sd.columns); h[0] = "id"; sd.columns = h

    d1 = pd.merge(bd, sd, how="left",  left_on="id", right_on="id") ## default: how="inner"
    # d1 = d1[d1[["length", "sequence"]].apply(lambda x: x[0] == len(x[1]), axis=1)]
    d1["read_sequence"] = d1["sequence"]
    d1["read_length"] = d1["read_sequence"].apply(lambda x: len(x)) 
    d1["sequence"] = d1[["sequence", "qstart", "qend"]].\
        apply(lambda x: x[0][(x[1]-1):x[2]], axis=1)

    if d1.shape[0] == 0:
        print("now rows for", s)
        continue

    d1["ss"] = d1[["sstart", "send"]].apply(lambda x: min(x[0], x[1]), axis=1)
    d1["se"] = d1[["sstart", "send"]].apply(lambda x: max(x[0], x[1]), axis=1)
    d1["strand"] = d1[["sstart", "ss"]].apply(lambda x: "+" if x[0] == x[1] else "-", axis=1)
    d1["sstart"], d1["send"] = d1["ss"], d1["se"]
    d1 = d1.drop(["ss", "se"], axis=1)
# d1["alin_seq"] = d1[["sequence", "strand"]].apply(lambda x:
#    x[0] if x[1]=="+" else str(Seq(x[0]).reverse_complement()), axis=1)

    h = ["id", "sequence", "strand", "sseqid", "pident", "mismatch", "gapopen", \
      "length", "sstart", "send", "evalue", "bitscore", "copy", "read_sequence", "read_length"]

    d2 = d1[h]
    h[3], h[8], h[9] = "ref", "start", "end"
    d2.columns = h
    d2 = d2.sort_values(["start"])
    d2.index.name = "id"
    d2.to_csv(out, sep="\t", index=False)
    print("    saved", out)
