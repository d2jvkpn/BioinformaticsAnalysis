# https://github.com/d2jvkpn/BioinformaticsAnalysis

import os
from Bio import SeqIO
import textwrap

for r in SeqIO.parse(os.sys.argv[1], "fasta"):
    print ('>' + r.id)
    codon = textwrap.wrap(str(r.seq[len(r)%3:]), 3)
    print( '\n'.join(codon))
