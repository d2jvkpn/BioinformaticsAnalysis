#! python3

import os, re, gzip

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']:
   print('''Convert KEGG pathway keg file to tsv format, usage:
    python3 keg2tsv.py  [.keg(.gz)]  [output.tsv]

project: https://github.com/d2jvkpn/BioinformaticsAnalysis
''')
   os.sys.exit(0)

keg, tsv = os.sys.argv[1:3]

if keg.endswith('.gz'):
    with gzip.open(keg, 'r') as f:
        kegtext = [ i.decode('utf8') for i in f.readlines()]
else:
    with open(keg, 'r') as f: kegtext = f.readlines()

TSV = open(tsv, 'w')

TSV.write ('\t'.join(['C_id', 'C_entry', 'C_name', 'gene_id', 'gene_name', \
'gene_description', 'A_id', 'A_name', 'B_id', 'B_name', 'KO', 'EC']) + '\n')

for line in kegtext:
    line = line.strip('\n')

    if len(line) in [0, 1]: continue
    if line[0] not in ['A', 'B', 'C', 'D']: continue 

    if line[0] == 'A':
        l1, L1 = re.split('\s+', line, 1)
    elif line[0] == 'B':
        l2, L2 = re.split('\s+', line, 2)[1:]
    elif line[0] == 'C' :
        l3, _ = re.split('\s+', line, 2)[1:]
        L3, k3 = _.rsplit(' ', 1) if '[' in _ else [_, '']
        TSV.write('\t'.join(['No.' + l3, k3.strip('[').strip(']'), L3, '', \
        '', '', l1, L1, 'No.' + l2, L2, '', '']) + '\n')
    else:
        c1, c2 = line.split('\t') if '\t' in line else [line, '; ']
        gi, _ = re.split('\s+', c1, 2)[1:]
        gn, gd = _.split('; ', 1) if '; ' in _ else ['', _]
        ko, ec = c2.split('; ', 1) if '; ' in c2 else [c2, '']
        
        TSV.write('\t'.join(['No.' + l3, '-', '-', gi, gn, gd, l1, '-', 
        'No.' + l2, '-', ko, ec]) + '\n')
