import os
## from glob import glob
import pandas as pd

# fs = glob("blastn/*.align.tsv") + glob("clean/*.tsv")
fs = os.sys.argv[1:]
for f in fs:
    if f.endswith(".length.tsv"): continue
    d1 = pd.read_csv(f, sep="\t")[["length", "copy"]]
    s1 = d1.groupby(["length"])["copy"].sum()
    s2 = pd.DataFrame({"count":s1.values}, index=s1.index)
    s2.to_csv(f.replace(".tsv", ".length.tsv"), sep="\t")