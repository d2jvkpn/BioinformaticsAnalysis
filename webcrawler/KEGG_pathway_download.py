import os, requests, time, subprocess, argparse, gzip
from multiprocessing import Pool

## get organisms code
# http://www.kegg.jp/kegg/catalog/org_list.html

__author__ = 'd2jvkpn'
__release__ = '2018-06-02'
__version__ = '0.1'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__lisence__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

#### 1. argument parser
parser = argparse.ArgumentParser ()

parser.add_argument ("org_codes", nargs = '+', 
help = 'organisms codes, E.g. hsa mmu ath')

parser.add_argument ("-p", dest="ParallelNum", type = int, default = 8,
help = "parallel request number, default: 8")

parser.add_argument ("-d", dest = "directory", default = "./",
help = "output directory, default: ./")

parser.add_argument ("-c", dest = "compress", type = bool, default = False,
help = "whether compress (gzip) keg file, default: False")

if len (os.sys.argv) == 1 or os.sys.argv[1] in ['-h', '--help']:
    print('Download KEGG pathway (.keg) file.\n')

    parser.print_help ()

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlisence: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __lisence__]
    print (_.format (*__))

    os.sys.exit (0)

args = parser.parse_args()


#### 2. define commands
def get_keg (org_code, compress, directory):
    keg_url = 'http://www.kegg.jp/kegg-bin/download_htext?htext=' + \
    org_code + '00001.keg&format=htext&filedir='

    resp = requests.get (keg_url)

    os.system ('mkdir -p ' + directory)

    fn = '{}/{}00001.keg'.format (directory, org_code)

    if compress == True:
        with gzip.open (fn + '.gz', 'wb') as f:
            f.write ( bytes (resp.text, 'utf8') )
    else:
        with open (fn , 'w') as f: f.write ( resp.text )


#### 3. initialize
ParallelNum = args.ParallelNum
if len (args.org_codes) < ParallelNum: ParallelNum = len(args.org_codes)


#### 4. run command
p = Pool (ParallelNum)

results = [p.apply_async (get_keg, \
args = (i, args.compress, args.directory)) for i in args.org_codes]

p.close ()
p.join ()
