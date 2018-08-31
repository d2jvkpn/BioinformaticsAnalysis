#! python3

__author__ = 'd2jvkpn'
__version__ = '1.7'
__release__ = '2018-08-27'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import os, gzip, sys
from collections import defaultdict
import pandas as pd

HELP = '''Convert ncbi gff3 format to ensembl gtf. Usage:
python3 ncbi_gff3_to_ensembl_gtf.py <input.gff.gz> <outputPrefix>'''

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']: 
    print (HELP)

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))
    os.sys.exit(0)

gff3, prefix = os.sys.argv[1:3]

####
GFF3 = gzip.open(gff3, 'rb')
records, protein, product  = [], defaultdict (str), defaultdict (str)

for line in GFF3:
    line = line.decode('utf8').strip()
    fd = line.split('\t')

    if fd[0].startswith('#') or len(fd)!=9: continue

    f9, d = fd[8].split(';'), defaultdict(str)

    for i in f9: ii = i.split('=', 1); d[ii[0]] = ii[1]

    if fd[2] == 'CDS' and 'protein_id' in d:
        protein[d['Parent']] = d['protein_id']

    if fd[2] in  ['exon', 'CDS'] and 'product' in d:
        product[d['Parent']] = d['product']

    if fd[2] in ['CDS', 'exon'] or int(fd[3]) >= int(fd[4]): continue

    if 'Parent' not in d and d['gbkey'].count('RNA') > 0:
        print ('Warning: invalid transcript record\t%s' % line, file=sys.stderr)
        continue

    if 'Parent' not in d and d['gbkey'] != 'Gene': continue

    if 'Parent' in d:
        _ = d['Name'] if d['Name'] != '' else d['transcript_id'] 
    else:
        _ = d['Name'] if d['Name'] != '' else d['gene']

    records.append ([d['ID'], fd[2], d['gbkey'], d['Parent'], _, \
    d['gene_biotype'], d['description'], d['Dbxref']])

GFF3.close()

####
M = pd.DataFrame.from_records (records, columns = ['id', 'feature', 'gbkey', \
'parent', 'name', 'gene_biotype', 'gene_description', 'Dbxref'])

M['product'] = [ product[i] for i in  M.iloc[:, 0]]
M['protein_id'] = [ protein[i] for i in  M.iloc[:, 0]]

M.to_csv (prefix + '.genes.tsv.gz', sep='\t', encoding='utf-8', 
index=False, compression='gzip')

# del (GFF3, records, protein, product, d, M); gc.collect()

M = pd.read_csv(prefix + '.genes.tsv.gz', sep='\t', index_col=0, 
usecols=[0, 3, 4, 5, 6, 8, 9])

M = M[[not i for i in M.index.duplicated()]]
M.fillna ('', inplace=True)

####
def toGtf(d):
    kv = []

    if 'gene_id' in d: kv.append('gene_id \"' + d['gene_id'] + '";')

    if 'transcript_id' in d:
        kv.append('transcript_id "' + d['transcript_id'] + '";')

    for k in d:
        if k in ['gene_id', 'transcript_id']: continue
        kv.append('%s "%s";' % (k, d[k]))

    return (' '.join(kv))

####
def Parser (d):
    del(d['gbkey'])

    if fd[2] == 'gene':
        d['gene_id'] = d['ID']
        d['gene_name'] = M.loc[d['gene_id'], 'name']
    elif fd[2] == 'transcript':
        d['transcript_id'] = d['ID']
        d['transcript_name'] = M.loc[d['ID'], 'name']
        d['gene_id'] = d['Parent']

        if not d['gene_id'].startswith('gene'):
        # primary transcript
            d['parent'] = d['Parent']
            d['gene_id'] = M.loc[d['gene_id'], 'parent']

        d['gene_name'] = M.loc[d['gene_id'], 'name']
        d['protein_id'] = M.loc[d['transcript_id'], 'protein_id']
    else:
        d['transcript_id'] = d['Parent']

        if d['transcript_id'].startswith('gene'):
        # not transcribed pseudogene
            d['gene_id'] = d['transcript_id'] 
        else:
            d['gene_id'] = M.loc[d['transcript_id'], 'parent']

        if not d['gene_id'].startswith('gene'):
        # primary transcript, V_gene_segment, C_gene_segment, J_gene_segment...
            d['parent'] = d['Parent']
            d['gene_id'] = M.loc[d['gene_id'], 'parent']

        d['gene_name'] = M.loc[d['gene_id'], 'name']

        if not d['transcript_id'].startswith('rna'):
        # V_gene_segment, C_gene_segment, J_gene_segment
            d['transcript_id'] = d['gene_id']

        d['transcript_name'] = M.loc[d['transcript_id'], 'name']

        if 'Dbxref' in d: del (d['Dbxref'])

GFF3 = gzip.open (gff3, 'rb')
GTF = gzip.open (prefix + '.gtf.gz', 'wb')

Attributions = ['ID', 'Parent', 'gbkey', 'gene_biotype', 'description', \
'product', 'protein_id', 'Dbxref']

for line in GFF3:
    fd = line.decode('utf8').strip().split('\t')

    if len(fd) != 9:
        print ('Warning: invalid record\t%s' % line.decode('utf8').strip(), file=sys.stderr)
        continue

    if int(fd[3]) > int(fd[4]): fd[3], fd[4] = fd[4], fd[3]

    f9, d = fd[8].split(';'), defaultdict(str)

    for _ in f9:
        i = _.split('=', 1)
        if i[0] in Attributions: d[i[0]] = i[1]

    if d['gbkey'] == 'Gene': fd[2] = 'gene'

    if d['gbkey'].count('RNA') > 0:
        d['transcript_biotype'] = d['gbkey']
        if fd[2] != 'exon': d['transcript_type'], fd[2] = fd[2], 'transcript'

    if fd[2] not in ['CDS', 'exon', 'transcript', 'gene']: continue

    try:
        Parser(d)
    except:
        print('Warning: "%s" missing attribution(s)' % d['ID'], file=sys.stderr)
        d['transcript_id'] = d['Parent'] if fd[2] == "exon" else d['ID']

    del (d['ID'])
    if 'Parent' in d: del (d['Parent'])

    for i in list(d.keys()):
        if d[i] == '': del(d[i])

    fd[8] = toGtf (d)
    GTF.write (bytes ('\t'.join(fd) + '\n', 'utf8'))

GFF3.close(); GTF.close()
