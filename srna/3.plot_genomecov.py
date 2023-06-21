import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

align, fai = os.sys.argv[1:3]

d0 = pd.read_csv(align, sep="\t")
xlen = pd.read_csv(fai, sep="\t", header=None, index_col=0).iloc[:, 0]
xlen = xlen.to_dict()
name = d0["ref"][0]
size = xlen[name]
prefix = align.replace(".tsv", "")

os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)


dp = d0.loc[d0["strand"] == "+", :].groupby(["start"])["copy"].sum()
dn = d0.loc[d0["strand"] == "-", :].groupby(["end"])["copy"].sum()

# dp = d0.loc[d0["strand"] == "+", ["sstart", "copy"]].copy()
# dp.index = dp["sstart"]
# dp = dp["copy"]
# dn = d0.loc[d0["strand"] == "-", ["sstart", "copy"]].copy()
# dn.index = dn["sstart"]
# dn = dn["copy"]


plt.bar(dp.index, dp.values, width=20, color='b')
plt.bar(dn.index, -dn.values, width=20, color='r')
plt.xlim(0, size)
#plt.ylim(-max(dn), dp)

s = 10**round(np.log10(size/10))
if size/s <=7: s = s/2

ticks=np.arange(0, size, s)
if len(ticks) > 15:
    plt.xticks(ticks, rotation=30)
else:
    plt.xticks(ticks)

r = max(dp) + max(dn)
s = 10**round(np.log10(r/10))
print(">>>", s, r)
if r/s < 7:
    s = s/2
elif r/s > 18:
    s = s*2

k = np.arange(0, max(dn), s)
m = (-k[::-1], np.arange(0, max(dp), s))

plt.yticks(np.concatenate(m, axis=None))

# this locator puts ticks at regular intervals
# loc = plticker.MultipleLocator(base=5)
# plt.xticks.set_major_locator(loc)

plt.title(name)
plt.xlabel("position")
plt.ylabel("reads count")

pdf = prefix + ".genomecov.pdf"
plt.savefig(pdf, bbox_inches = 'tight')
plt.close()
print("saved", pdf)
