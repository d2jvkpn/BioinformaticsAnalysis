#! /usr/bin/python3

__author__ = 'd2jvkpn'
__version__ = '0.2'
__release__ = '2018-07-28'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import requests, os, string, time
from bs4 import BeautifulSoup
from collections import OrderedDict

if len(os.sys.argv) == 1 or os.sys.argv[1] in ['-h', '--help']:
    print('''
Search species genome in NCBI by enter scientific name, get genome ftp
download by provide genome address.

Usage:
    python3 NCBIgenome.py search "species scientific name"
    e.g. "Mus musculus", Mus_musculus

    python3 NCBIgenome.py getftp "NCBI genome address"
    e.g. https://www.ncbi.nlm.nih.gov/genome/?term=Canis+lupus+familiaris[orgn]
''')


####
def formatSpeciesName(s):
    wds = s.replace('_', ' ').split()

    for i in range(len(wds)):
        a = False not in [i in string.ascii_letters for i in wds[i]]
        b = False not in [i in string.ascii_uppercase for i in wds[i]]
        if a and not b: wds[i] = wds[i].lower()

    wds[0] = wds[0].capitalize()
    return(' '.join(wds))


def search(name):
    _ = '+'.join (name.split()) + '%5Borgn%5D'
    link = 'https://www.ncbi.nlm.nih.gov/genome/?term=' + _

    query = requests.get (link)
    msg = 'NotFound' if 'No items found.' in query.text else link

    print(msg)
   
    return(1 if msg == 'NotFound' else 0)


def getftp(url):
    at = time.strftime ('%Y-%m-%d %H:%M:%S %z')
    r = requests.get(url)
    html_soup = BeautifulSoup(r.text, 'html.parser')
    Lineage = html_soup.find('span', class_='GenomeLineage').find_all('a')
    ln = [i.text  for i in Lineage[0::2]]
    txid = [i.get('href').split('/')[-1] for i in Lineage[0::2]]
    
    RGS = html_soup.find('div', class_='refgenome_sensor')
    _ = RGS.find_all ('span', class_='shifted')
    anchors = sum([i.find_all('a') for i in _], [])
    
    ftp = OrderedDict()
    
    for i in anchors:
        href = i.get('href')
        if href.startswith('ftp://') and os.path.basename(href) != '':
            ftp[i.text] = href
    
    loca = 'NCBI__' + ln[-1].replace(' ', '_') + '__'
    loca +=  list(ftp.values())[0].split('/')[-1].rsplit('_', 1)[0]

    os.system('mkdir -p ' + loca)
    F = open(loca + '/genome.infor.txt', 'w')
    SH = open(loca + '/download.sh', 'w')
    
    F.write('URL: %s\n' % url)
    F.write('Access time: %s\n' % at)
    
    if len(RGS.select('b')) == 3:
        F.write('Reference genome: ' + RGS.select('a')[0].text + '\n\n')
    else:
        F.write('\n')
    
    F.write('Lineage name:\n    %s\n\n' % ', '.join(ln))
    F.write('Lineage txid:\n    %s\n\n' % ', '.join(txid))
    
    prefix = list(ftp.values())[0].rsplit('/', 1)[0]
    SH.write('#! /bin/bash\n\n## URL: %s\n' % url)
    SH.write('## Species: %s\n## Access time: %s\n\n' % (ln[-1], at))
    SH.write('prefix=%s\n\n' % prefix)
    k=0
    wget = 'wget -c -O $(dirname $0)/{0} -c $prefix/{1} -o {0}.download.log &\n'
    
    for i in ftp:
        k += 1
        F.write('%s:\n    %s\n\n' % (i, ftp[i]))
        _ = ftp[i].replace(prefix + '/', '')
        SH.write(wget.format (_.split('_')[-1], _))
        SH.write('p%s=$! \n\n' % k)
    
    SH.write('wait ' + ' '.join(['$p' + str(i+1) for i in range(k)]) + '\n')
    F.close(); SH.close()
    
    print(loca)

# arg1 = Glycine_max
# arg1 = "Glycine max"
# arg1 = 'https://www.ncbi.nlm.nih.gov/genome/?term=Glycine+max%5Borgn%5D'

cmd, arg1 = os.sys.argv[1:3]

if cmd == 'search':
    os.sys.exit (search (formatSpeciesName(arg1)))
elif cmd == 'getftp':
    getftp (arg1)
else:
    pass
