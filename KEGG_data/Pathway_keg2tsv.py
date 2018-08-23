#! python3

import os, re, gzip

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']:
   print('''Convert KEGG pathway keg file to tsv format, usage:
    python3 keg2tsv.py  <.keg(.gz)>  <outputPrefix>

project: https://github.com/d2jvkpn/BioinformaticsAnalysis
''')
   os.sys.exit(0)

keg, prefix = os.sys.argv[1:3]

if keg.endswith('.gz'):
    with gzip.open(keg, 'r') as f:
        kegtext = [ i.decode('utf8') for i in f.readlines()]
else:
    with open(keg, 'r') as f: kegtext = f.readlines()

TSV = open(prefix + '.tsv', 'w')

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

TSV.close()


CMD = '''
awk 'BEGIN{FS=OFS="\t"} $2~"^PATH:"{ a[$1] = $2"\t"$3"\t"$8"\t"$10 }
a[$1] && $2!~"^PATH:"{split($11,x," "); KO=x[1];
if($12~"[[]EC:") {split($12,y," "); EC=y[length(y)]} else {EC=""};
print a[$1], $4, KO, EC}' $prefix.tsv |
awk 'BEGIN{FS=OFS="\t"} {a[$1]=$2"\t"$3"\t"$4; if(++x[$1"\t"$5]==1) g[$1]++; 
if(++x[$1"\t"$6]==1) k[$1]++; if($7!="" && ++x[$1"\t"$7]==1) e[$1]++}
END{for(i in a) {if (e[i]=="") e[i]=0; print i,a[i], g[i], k[i], e[i]} }' |
sort | awk 'BEGIN{FS=OFS="\t"; 
print "pathway", "pathway_name", "L1", "L2", "gene_count", "KO_count", "EC_count"}
{print}' > $prefix.summary.tsv

awk 'BEGIN{FS=OFS="\t"} $2~"^PATH:"{ a[$1] = $2"\t"$3"\t"$8"\t"$10 }
a[$1] && $2!~"^PATH:"{split($11,x," "); KO=x[1];
if($12~"[[]EC:") {split($12,y," "); EC=y[length(y)]} else {EC=""};
print a[$1], $4, KO, EC}' $prefix.tsv |
awk 'BEGIN{FS=OFS="\t"} {if(!idx[$4]) {n++; idx[$4]=n}
a[$4]=$3; if(++x[$4"\t"$5]==1) g[$4]++;
if(++x[$4"\t"$1]==1) p[$4]++;
if(++x[$4"\t"$6]==1) k[$4]++; if($7!="" && ++x[$4"\t"$7]==1) e[$4]++}
END{for(i in a) {if (e[i]=="") e[i]=0; print idx[i],i, a[i], p[i], g[i], k[i], e[i]} }' |
sort -n | awk 'BEGIN{FS=OFS="\t"; 
print "pathway_L2", "pathway_L1", "pathway_count", "gene_count", "KO_count", "EC_count"}
{print $2,$3,$4,$5,$6,$7}' > $prefix.classfication.tsv
'''

os.system ('prefix=%s\n' % prefix + CMD)
