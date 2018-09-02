#! python3

import os

def formatSpeciesName (s):
    import string

    wds = s.split()

    for i in range(len(wds)):
        a = False not in [i in string.ascii_letters for i in wds[i]]
        b = False not in [i in string.ascii_uppercase for i in wds[i]]
        if a and not b: wds[i] = wds[i].lower()

    wds[0] = wds[0].capitalize()
    return(' '.join (wds))

species = formatSpeciesName (os.sys.argv[1])
tsv = os.path.dirname (__file__) + '/data/KEGG_organism.tsv'
TSV = open (tsv, 'r')
ks = TSV.readline().strip('\n').split("\t")
Found = False
print ()

for line in TSV:
    fds = line.strip ().split("\t")
    if fds[2] == species:        
        for i in range (len (ks)): print (ks[i] + ": " + fds[i])
        Found = True; break

TSV.close ()
if not Found: print ("Not found")
print ()
