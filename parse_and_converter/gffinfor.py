#! /usr/bin/env python3

__author__ = 'd2jvkpn'
__version__ = '0.5'
__release__ = '2018-05-27'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import os, gzip
from collections import defaultdict

if len(os.sys.argv) != 4 or os.sys.argv[1] in ['-h', '--help']: 
    print ('Extract attribution(s) of feature(s) from gtf or gff3, usage:')
    _ = '<feature1,feature2...>  <attribution1,attribution2...>'
    print ('    python3 gffinfor.py  <input.gtf/input.gff3> ' + _  )
    print('\nOutput: stdout, tsv format.')

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))
    os.sys.exit(0)

gtf = os.sys.argv[1]
fea = os.sys.argv[2].split(',')
att = os.sys.argv[3].split(',')

GTF = gzip.open(gtf, 'rb') if gtf.endswith('gz') else open(gtf, 'r')
print ('\t'.join (att))

for _ in GTF:
    if gtf.endswith('gz'):
        fd = _.decode('utf8').strip('\n').strip(';').split('\t')
    else:
        fd = _.strip('\n').strip(';').split('\t')

    if fd[0].startswith('#') or len(fd) != 9 or fd[2] not in fea: continue

    d = defaultdict (str)
    delimiter, sep = (' ', '; ')  if fd[8].count('; ') > 0 else ('=', ';')

    for _ in fd[8].split(sep):
        i = _.split(delimiter, 1); d[i[0]] = i[1].strip('"')

    print (*[d[i] for i in att], sep='\t')

GTF.close()
