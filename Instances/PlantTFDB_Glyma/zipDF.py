import pandas as pd

Df = pd.read_csv("GLYMA0.tf.tsv", sep="\t", header=0)

md = Df.groupby(by="GeneID", axis=0)[Df.columns[1]].apply(lambda x: '; '.join(x))
md = md.to_frame()

for i in Df.columns[2:]:
    d = Df.groupby(by="GeneID", axis=0)[i].apply(lambda x: '; '.join(x))
    md[i] = d[md.index]

md.to_csv("GLYMA.tf.tsv", sep="\t")
