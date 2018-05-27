#! /usr/bin/env python3

__author__ = 'd2jvkpn'
__version__ = '0.8'
__release__ = '2018-05-27'
__project__ = 'https://github.com/d2jvkpn/GenomicProcess'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import os, gzip
from collections import defaultdict
import pandas as pd

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']: 
    print ('Convert ncbi gff3 format to ensembl gtf. Usage:')
    print ('python3 ncbi_gff3_to_ensembl_gtf.py <input.gff.gz> <outputPrefix>')

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))
    os.sys.exit(0)

gff3 = os.sys.argv[1]
tsv = os.sys.argv[2] + '.tsv'
gtf = os.sys.argv[2] + '.gtf.gz'

####
def toGtf(d):
    c9 = 'gene_id \"' + d['gene_id'] + '";'
    if 'transcript_id' in d:
       c9 += ' transcript_id "' + d['transcript_id'] + '";'

    for k in d:
        if k in ['gene_id', 'transcript_id']: continue
        c9 += ' %s "%s";' % (k, d[k])

    return (c9)


####
GFF3 = gzip.open(gff3, 'rb')

TSV = open(tsv, 'w')

TSV.write('\t'.join(['id', 'feature','gbkey','parent','name','gene_biotype', \
'description', 'product']) + '\n')

for _ in GFF3:
    fd = _.decode('utf8').strip().split('\t')

    if fd[0].startswith('#') or len(fd)!=9: continue

    f9=fd[8].split(';')
    d = defaultdict(str)

    for i in f9: ii = i.split('=', 1); d[ii[0]] = ii[1]

    if fd[2] in ['CDS', 'exon'] or int(fd[3]) >= int(fd[4]): continue
    if 'Parent' not in d and d['gbkey'] != 'Gene': continue

    if 'Parent' in d:
        _ = d['Name'] if d['Name'] != '' else d['transcript_id'] 
    else:
        _ = d['Name'] if d['Name'] != '' else d['gene']

    TSV.write('\t'.join([d['ID'], fd[2], d['gbkey'], d['Parent'], _, \
    d['gene_biotype'], d['description'], d['product']]) + '\n')

GFF3.close()
TSV.close()

####
M = pd.read_csv(tsv, sep='\t', index_col=0, usecols=[0,3,4,5])
M = M[[not i for i in M.index.duplicated()]]
M.fillna ('', inplace=True)

GFF3 = gzip.open(gff3, 'rb')
GTF = gzip.open(gtf, 'wb')

Attributions = ['ID', 'Parent', 'protein_id', 'gbkey', 'gene_biotype',
'description', 'product']

for _ in GFF3:
    fd = _.decode('utf8').strip().split('\t')

    if fd[0].startswith('#') or len(fd)!=9: continue

    f9=fd[8].split(';')
    d = {}

    for _ in f9:
        i = _.split('=', 1)
        if i[0] in Attributions: d[i[0]] = i[1]

    if int(fd[3]) >= int(fd[4]) or 'gbkey' not in d: continue

    if d['gbkey'] == 'Gene': fd[2] = 'gene'

    if d['gbkey'].count('RNA') > 0 and fd[2] != 'exon':
        d['transcript_type'] = fd[2]
        d['transcript_biotype'] = d['gbkey']
        fd[2] = 'transcript'

    if fd[2] not in ['CDS', 'exon', 'transcript', 'gene']: continue

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

        if 'product' in d: del(d['product'])


    del(d['ID'])
    if 'Parent' in d: del(d['Parent'])

    for i in list(d.keys()):
        if d[i] == '': del(d[i])

    fd[8] = toGtf (d)
    GTF.write (bytes ('\t'.join(fd) + '\n', 'utf8'))

GFF3.close()
GTF.close()
