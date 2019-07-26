import os
from Bio import SeqIO

__author__ = 'd2jvkpn'
__version__ = '0.2'
__release__ = '2019-07-26'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__lisense__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

if len(os.sys.argv) != 4:
    print("Split fasta into parts, usage:")
    print("  $ python3 fasta-splitter.py  <input.fasta>  <n>  <output_prefix>")
    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisense__]
    print (_.format (*__))
    os.sys.exit(2)

fa, n, prefix = os.sys.argv[1:4]
n = int(n)
if n < 1: print("invalid parts number:", n); os.sys.exit(1)
# fa, n, prefix = "mRNA.fa", 20, "out"

m = 0
sfa = SeqIO.parse(fa, "fasta") 
for r in sfa: m+=1
sfa.close()

n = min(n, m)
q, r = divmod (m, n)
qs = [q]*n
for i in range(r): qs[i] +=1
ns = [ prefix + ".part-%d.fasta" % (i+1) for i in range(n)]
os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)

i, x = 0, 0
sfa = SeqIO.parse(fa, "fasta")
for r in sfa:
    if i == 0: f = open(ns[x], "w")
    SeqIO.write(r, f, 'fasta'); i += 1
    if i == qs[x]:
         f.close()
         print("saved {}, {} records".format(ns[x], qs[x]))
         i, x = 0, x + 1

sfa.close()
