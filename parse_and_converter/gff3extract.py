import os, gzip, argparse


__author__ = 'd2jvkpn'
__version__ = '0.9'
__release__ = '2018-06-20'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

parser = argparse.ArgumentParser(description='Extract transcript records')

parser.add_argument("gff3", help="input gff3 file")

parser.add_argument("-F", dest="Feature", nargs='+', required=True, 
help="transcript type(s) and attribution(s), " +
"E.g \"-F mRNA ncRNA:ncrna_class=lncRNA\"")

parser.add_argument("-A", dest="Attri", nargs='+', default=None, 
help="output attribution(s) of selected transcript(s)")

parser.add_argument("-a", dest="attri", nargs='+', default=None, 
help="output attribution(s) of exon")

if len(os.sys.argv) == 1:
    parser.print_help()
    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))
    os.sys.exit(0)

args = parser.parse_args()
MFeature = {}

for i in args.Feature:
    ii = i.split(":")

    try:
        MFeature[ ii[0] ] = ii[1].split(",")
    except:
        MFeature[ ii[0] ] = None

def extract ( Fx9 ):
    global k; f9 = Fx9[8].split(';'); t = Fx9[2]

    if t in MFeature:
        k = 1
        if MFeature[t] != None and not any(x in MFeature[t] for x in f9):
            k=0

    if t != 'exon' and t not in MFeature: k = 0
    if k != 1: return 0

    atr = args.Attri if t in MFeature else args.attri
    F9 = ""

    if atr != None:
        for i in f9:
            if i.split(sep='=', maxsplit=1)[0] in atr: F9 += i + ';'

        Fx9[8] = F9.strip(';')

    print('\t'.join(Fx9))


if args.gff3.endswith('.gz'):
    f = gzip.open(args.gff3, 'r')
    for line in f:
        line = line.decode("utf8").rstrip()
        if not line.startswith("#"): extract( line.split('\t') )
else:
    f = open(args.gff3, 'r')
    for line in f:
        line = line.rstrip()
        if not line.startswith("#"): extract( line.split('\t') )

f.close()
