import os
from Bio import SeqIO

__author__ = 'd2jvkpn'
__version__ = '0.3'
__release__ = '2019-08-31'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__lisense__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

HELP = '''Split fasta into parts, usage:
  $ python3 fasta-splitter.py  <input.fasta>  <n>  <output_prefix>

author: {}
version: {}
release: {}
project: {}
lisence: {}
'''

if len(os.sys.argv) != 4:
    _ = [__author__,  __version__, __release__, __project__, __lisense__]
    print (HELP.format (*_), file=os.sys.stderr)
    os.sys.exit(2)

fa, n, prefix = os.sys.argv[1:4]
n = int(n)
if n < 1: os.sys.exit("invalid parts number: %d" % n)
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
print("Split {} into {} parts...".format(fa, n), file=os.sys.stderr)
i, x, tmp, sfa = 0, 0, "saved {}, {} records", SeqIO.parse(fa, "fasta")
for r in sfa:
    if i == 0: f = open(ns[x], "w")
    SeqIO.write(r, f, 'fasta'); i += 1
    if i == qs[x]:
        f.close()
        print(tmp.format(ns[x], qs[x]), file=os.sys.stderr)
        i, x = 0, x + 1

sfa.close()

print(n)
