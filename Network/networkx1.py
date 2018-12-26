__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'

import os
import pandas as pd
import networkx as nx

itsv, prefix = os.sys.argv[1], os.sys.argv[2]

d = pd.read_csv(itsv, sep="\t", header = 0)

dedges = [(i[1], i[2]) for i in d.iloc[:, [0, 1]].to_records()]
# G = nx.from_pandas_edgelist(df=d, source="#node1", target="node2")

G = nx.Graph()
G.add_edges_from(dedges)

dx = pd.DataFrame.from_records(list(G.degree), columns=["gene", "degree"])

dt = nx.betweenness_centrality(G)
dx["betweenness"] = [dt[g] for g in dx.gene]

subg = sorted(nx.connected_components(G), key=len, reverse=True)

for i in range(len(subg)):
    tsv = prefix + ".subgrap_{}.tsv".format(i+1)
    d1 = d.loc[ [(j[0] in subg[i] and j[1] in subg[i]) for j in dedges], :]
    d1.to_csv(tsv, sep="\t", index=False)
    print("Saved", tsv)


def InSubgrpahs(g, subg):
    x = []
    for i in range(len(subg)):
        if g in subg[i]: x.append(i+1)

    return ', '.join([str(j) for j in x])
    
dx["subgraph"] = [InSubgrpahs(g, subg) for g in dx.gene]

tsv = prefix + ".network_infor.tsv"
dx.to_csv(tsv, sep="\t", index=False)
print("Saved", tsv)

# max(nx.connected_components(G),key=len)
# [len(c) for c in sorted(nx.weakly_connected_components(G), key=len, reverse=True)]
# nx.strongly_connected_components(G) 
