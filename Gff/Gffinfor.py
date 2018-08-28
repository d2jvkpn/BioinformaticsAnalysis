#! python3

__author__ = 'd2jvkpn'
__version__ = '0.6'
__release__ = '2018-08-20'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import os, gzip
from collections import defaultdict

HELP = '''
Summary sequences, sources, types of gff/gtf(.gz):
    python3 Gffinfor.py  <gff>

Summary types' attributions gff/gtf:
    python3 Gffinfor.py  <gff>  <type1,type2...>

Extract attributions and Dbxref (tsv format) from gff/gtf:
("type" can be set as "")
    python3 Gffinfor.py  <gff>  <type1,type2...> <attr1,attr2,attr3>  [dbxref1,dbxref2...]
'''

if len (os.sys.argv) == 1 or os.sys.argv [1] in ['-h', '--help']:
    print (HELP)

    _ = 'author: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))

    os.sys.exit (0)

gtf = os.sys.argv [1]

CountC, CountS, CountT = defaultdict (int), defaultdict (int), defaultdict (int)
TypeA, TypeV = defaultdict (int), defaultdict (set)

if len (os.sys.argv) >= 3: Type = os.sys.argv [2].split (',')
if len (os.sys.argv) >= 4: Attr = os.sys.argv [3].split (',')
if len (os.sys.argv) == 4: print (*Attr, sep="\t")

if len (os.sys.argv) >= 5:
    Dbxref = os.sys.argv [4].split (',')
    print (*(Attr + Dbxref), sep="\t")

GTF = gzip.open (gtf, 'rb') if gtf.endswith ('gz') else open (gtf, 'r')

if gtf.endswith('.gtf') or gtf.endswith('.gtf.gz'):
    delimiter, sep = ' ', '; '
else:
    delimiter, sep = '=', ';'

for _ in GTF:
    if gtf.endswith ('gz'):
        fd = _.decode ('utf8').strip ('\n').strip (';').split ('\t')
    else:
        fd = _.strip ('\n').strip (';').split ('\t')

    if fd [0].startswith ('#') or len (fd) != 9: continue

    if len (os.sys.argv) == 2:
        CountC [fd [0]] += 1
        CountS  ["source " + fd [1]] += 1
        CountT["type " + fd [2]] += 1
        continue

    d = defaultdict (str)
    for _ in fd [8].split (sep):
        i = _.split (delimiter, 1); d [i [0]] = i [1].strip ('"')

    if len (Type) > 0 and fd [2] not in Type: continue

    if len (os.sys.argv) == 3:
        for i in d: k = fd [2] + "\t" + i; TypeA [k] += 1; TypeV [k].add (d [i])
        continue

    if len (os.sys.argv) == 4:
        print (*[d [i] for i in Attr], sep='\t')
        continue

    dbref = defaultdict (str)

    if d ['Dbxref'] != "":
        for i in d ['Dbxref'].split (','):
            ii = i.split (':', 1); dbref [ii [0]] = ii [1]

    print (* ([d [i] for i in Attr] + [dbref [i] for i in Dbxref]), sep='\t')

GTF.close ()

FMT = "echo '%s' | column -t -s $'\t' | awk 'NR==1{print \"0    \"$0; print 1} \
NR>1{print \"2    \"$0}' | sort | cut --complement -c 1"

if len (os.sys.argv) == 2:
    sl = ["NAME\tCOUNT", "Seqences\t" + str (len (CountC))]
    for i in CountS: sl.append (i + "\t" + str (CountS [i]))
    for i in CountT: sl.append (i + "\t" + str (CountT [i]))

    os.system (FMT % '\n'.join (sl))
elif len (os.sys.argv) == 3:
    sl = ['TYPE\tATTRIBUTION\tTOTAL\tUNIQUE']

    for i in TypeA:
        sl.append (i + "\t" + str( TypeA [i]) + "\t" + str (len (TypeV [i])))

    os.system (FMT % '\n'.join (sl))
else:
    pass
