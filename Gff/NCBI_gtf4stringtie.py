import os, gzip, sys

if len(os.sys.argv) == 1 or os.sys.argv[1] in ['-h', 'help']:
    print ("Replace values of gene_id and transcript_id (NCBI gtf) with " + \
    "gene_name and transcript_name\nNote: output gtf will contains duplicate ids.")

    os.sys.exit(0)

gtf = os.sys.argv[1]
fr = gzip.open (gtf, 'rb')

for line in fr:
    fd = line.decode('utf8').strip(';\n').split('\t')

    if fd[0].startswith('#') or len(fd) != 9  or int(fd[3]) >=  int(fd[4]) or \
    fd[2] not in ['gene','transcript', 'exon', 'CDS']:
        continue

    d, kv = {}, []

    for _ in fd [8].split ('; '):
        i = _.split (' ', 1); d [i [0]] = i [1].strip ('"')

    if 'gene_name' in d: d['gene_id'] = d['gene_name']

    if 'transcript_name' in d:
        d['transcript_id'] = d['transcript_name']
        d.pop('transcript_name')

    a = 'gene_id' in d and 'transcript_id' in d

    if fd[2] in ['transcript', 'exon'] and not a:
        print ('Warning: skip record\t%s' % line.decode('utf8'), file=sys.stderr)
        continue

    if 'gene_id' in d: kv.append('gene_id "%s";' % d['gene_id'])

    if 'transcript_id' in d:
        kv.append ('transcript_id "%s";' % d['transcript_id'])

    for k in d:
        if k in ['gene_id', 'transcript_id']: continue
        kv.append('%s "%s";' % (k, d[k]))

    fd[8] = ' '.join(kv)
    print('\t'.join(fd))

fr.close()
