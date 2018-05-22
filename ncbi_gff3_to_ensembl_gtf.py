#! /usr/bin/env python3

__author__ = 'd2jvkpn'
__version__ = '0.4'
__release__ = '2018-05-20'
__project__ = 'https://github.com/d2jvkpn/GenomicProcess'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import os, gzip
from collections import defaultdict
import pandas as pd

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']: 
    print ('Convert ncbi gff3 format to ensembl gtf. Usage:')
    print ('   python3 ncbi_gff3_to_ensembl_gtf.py <input.gff3> <output.gtf>')
    print ('\nNote: file "genomic.transcription.tsv" will be created, "pandas" is required.')

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))
    os.sys.exit(0)

gff3 = os.sys.argv[1]
gtf = os.sys.argv[2]

####
def toGtf(d):
    c9 = 'gene_id \"' + d['gene_id'] + '";'
    if 'transcript_id' in d:
       c9 += ' transcript_id "' + d['transcript_id'] + '";'

    for k in d:
        if k in ['gene_id', 'transcript_id']: continue
        c9 += ' %s "%s";' % (k, d[k])

    return (c9)

def Mapper(d, M, t):
    Q = (d['ID'] if t == 'transcript' else d['Parent'])
    d['transcript_id'] = Q

    for _ in [1,2,3]:
        if not d['transcript_id'].startswith('rna'): d['transcript_id'] =  Q
        Q = M.ix[Q, 'parent']
        if Q == '': Q = d['transcript_id']; break
        if M.ix[Q, 'gbkey'] == 'Gene': d['gene_id'] = Q; break

    if _ == 2: d['parent'] = M.ix[d['transcript_id'], 'parent']

    if 'gene_id' not in d: d['gene_id'] = d['transcript_id']
    d['gene_biotype'] = M.ix[Q, 'gene_biotype']

    del(d['ID']); del(d['Parent'])


####

GFF3 = gzip.open(gff3, 'r') if gff3.endswith('gz') else open(gff3, 'r')
tsv = open('genomic.transcription.tsv', 'w')
tsv.write('id\tfeature\tgbkey\tparent\tname\tgene_biotype\n')

for _ in GFF3:
    fd = _.decode("utf8").replace('\n', '').split('\t')
    if fd[0].startswith('#') or len(fd)!=9: continue

    f9=fd[8].split(';')
    d = defaultdict(str)

    for i in f9: ii = i.split('=', 1); d[ii[0]] = ii[1]

    if fd[2] in ['CDS', 'exon'] or int(fd[3]) >= int(fd[4]): continue

    if 'Parent' in d or d['gbkey'] == 'Gene':    
        tsv.write ('\t'.join([d['ID'], fd[2], d['gbkey'], d['Parent'], \
        d['Name'], d['gene_biotype']]) + '\n')


tsv.close()
GFF3.close()


####
M = pd.read_csv('genomic.transcription.tsv', sep='\t', index_col=0)
M.fillna ('', inplace=True)

GFF3 = gzip.open(gff3, 'r') if gff3.endswith('gz') else open(gff3, 'r')
GTF = open(gtf, 'w')

for _ in GFF3:
    fd = _.decode("utf8").replace('\n', '').split('\t')
    if fd[0].startswith('#') or len(fd)!=9: continue

    f9=fd[8].split(';')
    d = {}
    for i in f9: ii = i.split('=', 1); d[ii[0]] = ii[1]

    if int(fd[3]) >= int(fd[4]) or 'gbkey' not in d: continue
    if 'Parent' not in d and d['gbkey'] != 'Gene': continue

    if d['gbkey'].count('RNA') > 0 and fd[2] != 'exon':
        d['feature'] = fd[2]; fd[2] = 'transcript'

    if d['gbkey'] == 'Gene':
        fd[2] = 'gene'; d['gene_id'] = d['ID']; del(d['ID'])

    if fd[2] not in ['CDS', 'exon', 'transcript', 'gene']: continue

    if 'gene' in d: d['gene_name'] = d['gene']; del(d['gene'])

    if 'transcript_id' in d:
        d['transcript_name'] = d['transcript_id']; del(d['transcript_id'])

    if fd[2] != 'gene': Mapper(d, M, fd[2])

    fd[8] = toGtf (d)

    GTF.write ('\t'.join(fd) + '\n')

GFF3.close()
GTF.close()
