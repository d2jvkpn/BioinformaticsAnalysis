#! /usr/bin/env python3

__author__ = 'd2jvkpn'
__version__ = '0.3'
__release__ = '2018-05-20'
__project__ = 'https://github.com/d2jvkpn/GenomicProcess'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'


from collections import defaultdict
import os, gzip

if len(os.sys.argv) != 3 or os.sys.argv[1] in ['-h', '--help']: 
    print ('Convert ncbi gff3 format to ensembl gtf. Usage:')
    print ('   python3 ncbi_gff3_to_ensembl_gtf.py <input.gff3> <output.gtf>')
    print ('\nNote: transcription mapping file "genomic.transcription.tsv" will be created.')

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))
    os.sys.exit(0)

gff3 = os.sys.argv[1]
gtf = os.sys.argv[2]

tsv = open('genomic.transcription.tsv', 'w')

f1 = gzip.open(gff3, 'r') if gff3.endswith('gz') else open(gff3, 'r')

tsv.write('Feature\tID\tParent\tName\tType\n')

for _ in f1:
    try:
       line = _.decode("utf8").strip()
    except:
       line = _.strip()

    if line.startswith('#'): continue
    fd=line.split('\t')
    if int(fd[3]) >= int(fd[4]):
        print('\nError: start position is greater than or equal to end position at:')
        print('%s' % line)
    f9=fd[8].split(';')
    d = defaultdict(str)
    for i in f9:
        ii=i.split('=', 1)
        d[ii[0]] = ii[1]

    if d['gbkey'].count('RNA') > 0 and fd[2] != 'exon':
        tsv.write('\t'.join([fd[2], d['ID'], d['Parent'], \
        d['transcript_id'], d['gbkey']]) + '\n')
        continue

    if d['gbkey'] == 'Gene':
        tsv.write('\t'.join(['Gene', d['ID'], '', d['Name'], d['gene_biotype']]) + '\n')


tsv.close()
f1.close()

tp=defaultdict(str)
nm=defaultdict(str)

with open('genomic.transcription.tsv', 'r') as _:
    for line in _:
        fd = line.replace('\n', '').split('\t')
#        if fd[0] in ['CDS', 'exon']: tl[fd[1]] = fd[2]
        if fd[0] != 'Gene':
            tp[fd[1]] = fd[1] if fd[2] == '' else fd[2]
            nm[fd[1]] = fd[3]
        if fd[0] == 'Gene': nm[fd[1]] = fd[3]


def gF9(d):
    ID = d['ID']
    attr = 'gene_id "' + ID + '"; '
    if nm[ID] != '' : attr += 'gene_name "' + nm[ID] +'"; ' 
    attr += 'gene_biotype "' + d['gene_biotype']  + '"; '

    if d['description'] != '':
        attr += 'description "' + d['description']  + '"; '

    if d['Dbxref'] != '':
        attr += 'Dbxref "' + d['Dbxref']  + '"; '

    return (attr)

def rF9(d):
    ID = d['ID']
    Parent = d['Parent']
    attr = 'transcript_id "' + ID + '"; '
    if Parent.startswith('rna'):
        attr += 'gene_id "' + tp[Parent] + '"; '
    else:
        attr += 'gene_id "' + Parent + '"; '

    if d['transcript_id'] != '' :
        attr += 'transcript_name "' + d['transcript_id'] + '"; '
    if nm[tp[Parent]] != '':
        attr += 'gene_name "' + nm[tp[d['ID']]] + '"; '
    if d['product'] != '':
        attr += 'product "' + d['product'] +'"; '

    return (attr)


def eF9(d):
    Parent = d['Parent']

    attr = 'transcript_id "' + Parent + '"; '

    if Parent.startswith('gene') == '':
        attr = 'gene_id "' + Parent + '"; '
    elif tp[tp[Parent]] != '':
        attr += 'gene_id "' + tp[tp[Parent]] + '"; '
        attr += 'primary_transcript "' + tp[Parent] + '"; '
    else:
        attr += 'gene_id "' + tp[Parent] + '"; '

    if d['transcript_id'] != '' :
        attr += 'transcript_name "' + d['transcript_id'] + '"; '

    if d['gene'] != '' :
        attr += 'gene_name "' + d['gene'] + '"; '

    return(attr)


f1 = gzip.open(gff3, 'r') if gff3.endswith('gz') else open(gff3, 'r')
f2 = open(gtf, 'w')

for _ in f1:
    try:
       line = _.decode("utf8").strip()
    except:
       line = _.strip()

    if line.startswith('#'): continue
    fd=line.replace('\n', '').split('\t')

    if int(fd[3]) >= int(fd[4]): continue

    f9=fd[8].split(';')
    d = defaultdict(str)
    for i in f9: ii=i.split('=', 1); d[ii[0]] = ii[1]

    if d['gbkey'] not in ['CDS', 'exon', 'Gene'] and \
    d['gbkey'].count('RNA') == 0: continue

    if d['gbkey'].count('RNA') > 0 and fd[2] != 'exon': fd[2] = 'transcript'

    fd[8] = 'gbkey "'+ d['gbkey']+ '";'

    if d['gbkey'] == 'Gene':
        fd[8] = gF9(d) + fd[8]
        f2.write('\t'.join(fd) + '\n')
        continue

    if d['gbkey'].count('RNA') > 0 and fd[2] != 'exon':
        fd[8] = rF9(d) + fd[8]
        f2.write('\t'.join(fd) + '\n')
        continue

    if d['gbkey'] == 'CDS' or fd[2] == 'exon':
        fd[8] = eF9(d) + fd[8]
        f2.write('\t'.join(fd) + '\n')

f1.close
f2.close()
