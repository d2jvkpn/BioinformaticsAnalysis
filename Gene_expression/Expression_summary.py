#! python3

import os
import pandas as pd
import numpy as np

tsv, prefix, threshold = os.sys.args[1], os.sys.args[2], float(os.sys.args[3])

expr = pd.read_csv(tsv, sep='\t', index_col=0)

####
top = []

for i in expr.columns:
    x = expr.loc[:, i].nlargest(10)
    top.append (list (x.index.values))
    top.append (x.tolist())

topdf = pd.DataFrame.from_records(top)

topdf.index = sum([[i + "__gene", i + "__value"] for i in expr.columns], [])
topdf.columns = [ "top_" + str(i) for i in range(1,11)]

topdf.to_csv (prefix + ".top10.tsv", sep="\t")

####
expr[expr < threshold] = pd.np.nan

d = expr.describe().T
d["median"] = expr.apply (lambda x: np.median(x[~x.isnull()]), axis=0)[d.index]
d = round(d, 2)

d.loc[:, "count"] = d.loc[:, "count"].astype(int)
d.columns = ["count", "mean", "std", "min", "Q1", "Q2", "Q3", "max", "median"]
d = d.loc[:, ["count", "mean", "median", "std", "min", "Q1", "Q2", "Q3", "max"]]

d.to_csv(prefix + ".summary.tsv", sep="\t")
