#! /bin/env python3

__author__ = 'd2jvkpn'
__version__ = '0.4'
__release__ = '2018-07-28'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

import os, requests, time, string
from urllib.parse import urlparse
from ftplib import FTP
import pandas as pd
from bs4 import BeautifulSoup
from biomart import BiomartServer

'''
http://asia.ensembl.org/index.html
http://metazoa.ensembl.org/index.html
http://plants.ensembl.org/index.html
http://fungi.ensembl.org/index.html
http://bacteria.ensembl.org/index.html
http://protists.ensembl.org/index.html

http://asia.ensembl.org/info/about/species.html
http://metazoa.ensembl.org/species.html
http://plants.ensembl.org/species.html
http://fungi.ensembl.org/species.html
http://bacteria.ensembl.org/species.html
http://protists.ensembl.org/species.html
'''

if len(os.sys.argv) == 1 or os.sys.argv[1] in ['-h', '--help']:
    print('''
Search species genome in Ensembl by enter scientific name, get genome ftp
download links and achive gene annotation (GO, kegg_enzyme, entrez) using 
biomart by provide Ensembl genome address.

Usage:
    python3 EnsemblGenome.py search "species scientific name"
    e.g. "Mus musculus", Mus_musculus

    python3 EnsemblGenome.py getftp/biomart "ensembl genome address"
    e.g. http://asia.ensembl.org/Mus_musculus/Info/Index

Note:
    Please use Python3.6 or higher.
''')

    _ = '\nauthor: {}\nversion: {}\nrelease: {}\nproject: {}\nlicense: {}\n'
    __ = [__author__,  __version__, __release__, __project__, __license__]
    print (_.format (*__))

    os.sys.exit(0)


####
def formatSpeciesName(s):
    wds = s.replace('_', ' ').split()

    for i in range(len(wds)):
        a = False not in [i in string.ascii_letters for i in wds[i]]
        b = False not in [i in string.ascii_uppercase for i in wds[i]]
        if a and not b: wds[i] = wds[i].lower()

    wds[0] = wds[0].capitalize()
    return(' '.join(wds))


def query(url):  
    query = requests.get(url)
    if not query.ok: os.sys.exit('Failed to request "%s"' % url)
    bs = BeautifulSoup(query.text, 'html.parser')

    _  = bs.find('span', class_ = 'header').select('a')[0]
    version = _.text.strip(')').split(' (')[1]

    dnaFTP = ''

    for i in bs.find_all('a', class_='nodeco'):
        if i.get('href').startswith('') \
        and i.text == 'Download DNA sequence': dnaFTP = i.get('href')

    netloc = urlparse(dnaFTP).netloc
    path = urlparse(dnaFTP).path

    _ = dnaFTP.replace('/dna/', '/')
    ensembl = _.split('/')[-4].replace('release', 'Ensembl')
    ScentificName = _.split('/')[-2]
    ScentificName = ScentificName[0].upper() + ScentificName[1:]

    loca = '__'.join([ensembl, ScentificName, version])
    return(netloc, path, loca)


def getftp(netloc, path, loca, url):
    at = time.strftime('%Y-%m-%d %H:%M:%S %z')

    ensembl, ScentificName, version = loca.split('__')

    ftp = FTP(netloc)
    ftp.login()
    
    for i in ftp.nlst(path):
        if i.endswith('.dna_sm.toplevel.fa.gz'): dna = 'ftp://' + netloc + i
    
    for i in ftp.nlst(path.replace('/dna/', '/cdna/')):
        if i.endswith('.all.fa.gz'): cdna = 'ftp://' + netloc + i
    
    for i in ftp.nlst(path.replace('/dna/', '/ncrna/')):
        if i.endswith('.ncrna.fa.gz'): ncrna = 'ftp://' + netloc + i
    
    for i in ftp.nlst(path.replace('/dna/', '/pep/')):
        if i.endswith('.pep.all.fa.gz'): pep = 'ftp://' + netloc + i
    
    for i in ftp.nlst(path.replace('/dna/', '').replace('/fasta/', '/gtf/')):
        if i.endswith('.gtf.gz') and not i.endswith('.abinitio.gtf.gz'):
            gtf = 'ftp://' + netloc + i
    
    ftp.close()
    
    os.system('mkdir -p %s'  % loca)
    
    with open (loca + '/genome.infor.txt', 'w') as f:
        f.write('URL: %s\n' % url)
        f.write('Acess time: %s\n' % at)
        f.write('Scentific name: %s\n' % ScentificName.replace('_', ' '))
        f.write('Assembly version: %s\n' % version)
        f.write('Ensembl version: %s\n\n' % ensembl)
        f.write('DNA fasta:\n    %s\n\n' % dna)
        f.write('cdna fasta:\n    %s\n\n' % cdna)
        f.write('ncrna fasta:\n    %s\n\n' % ncrna)
        f.write('pep fasta:\n    %s\n\n' % pep)
        f.write('annotation gtf:\n    %s\n' % gtf)
    
    wget = 'wget -c -O {0} {1} -o {0}.download.log &\np{2}=$!\n\n'
    
    with open (loca + '/download.sh', 'w') as f:
        f.write('#! /bin/bash\n\n## URL: %s\n' % url)
        f.write('## Species: %s\n' % ScentificName.replace('_', ' '))
        f.write('## Acess time: %s\n\n' % at)
        f.write(wget.format('genomic.fa.gz', dna, 1))
        f.write(wget.format('cdna.fa.gz', cdna, 2))
        f.write(wget.format('ncrna.fa.gz', ncrna, 3))
        f.write(wget.format('pep.fa.gz', pep, 4))
        f.write(wget.format('genomic.gtf.gz', gtf, 5))
        f.write('wait $p1 $p2 $p3 $p4 $p5\n')

    print(loca)


def biomart_anno(url, loca):
    urlp = urlparse(url)
    
    species = urlp.path.split('/')[1]
    code = species.split('_')[0][0].lower() + species.split('_')[1]
    
    server = BiomartServer('%s://%s/biomart' % (urlp.scheme, urlp.netloc))
    datasets = server.datasets

    print("Connected to Ensembl biomart...")
    
    _ = ['metazoa', 'plants', 'fungi', 'bacteria', 'protists']
    dn = code + ('_eg_gene'  if urlp.netloc.split('.')[0] in _ else '_gene_ensembl')
    ds = datasets[dn]

    #### ds.attribute_pages
    os.system('mkdir -p %s' % loca)
    
    s1 = ds.search({'attributes': ['ensembl_gene_id', 'go_id']})
    
    gene2go = pd.DataFrame.from_records(
      [str(i, encoding = 'utf-8').split('\t') for i in s1.iter_lines()], 
      columns = ['gene', 'GO_id'])
    
    gene2go = gene2go.loc[gene2go['GO_id'] != '', :]
    gene2go.drop_duplicates(inplace=True)
    gene2go.to_csv(loca + '/gene2go.tsv', sep='\t', index=False)

    ####
    s2 = ds.search({'attributes': ['ensembl_gene_id', 'entrezgene']})
    
    gene2entrez = pd.DataFrame.from_records(
      [str(i, encoding = 'utf-8').split('\t') for i in s2.iter_lines()], 
      columns = ['gene', 'entrez'])
    
    gene2entrez = gene2entrez.loc[gene2entrez['entrez'] != '', :]
    gene2entrez.drop_duplicates(inplace=True)
    gene2entrez.to_csv(loca + '/gene2entrez.tsv', sep='\t', index=False)

    ####
    s3 = ds.search({'attributes': ['ensembl_gene_id', 'kegg_enzyme']})
    
    gene2kegg = pd.DataFrame.from_records(
      [str(i, encoding = 'utf-8').split('\t') for i in s3.iter_lines()], 
      columns = ['gene', 'kegg_enzyme'])
    
    gene2kegg = gene2kegg.loc[gene2kegg['kegg_enzyme'] != '', :]
    gene2kegg.to_csv(loca + '/gene2kegg.tsv', sep='\t', index=False)

    ####
    s4 = ds.search({'attributes': ['ensembl_gene_id', 'gene_biotype', \
    'external_gene_name', 'description']})
    
    gene_infor = pd.DataFrame.from_records(
      [str(i, encoding = 'utf-8').split('\t') for i in s4.iter_lines()], 
      columns = ['gene', 'gene_biotype', 'gene_name', 'gene_description'])

    ####
    g = gene2go.groupby('gene')['GO_id'].apply(lambda x: ', '.join(x))
    k = gene2kegg.groupby('gene')['kegg_enzyme'].apply(lambda x: ', '.join(x))
    e = gene2entrez.groupby('gene')['entrez'].apply(lambda x: ', '.join(x))
    
    gene_infor['GO_id'] = [ g[i] if i in g else '' for i in gene_infor['gene']]
    gene_infor['kegg_enzyme'] = [ k[i] if i in k else '' for i in gene_infor['gene']]
    gene_infor['entrez'] = [ e[i] if i in e else '' for i in gene_infor['gene']]
    
    gene_infor.to_csv(loca + '/gene.infor.tsv', sep='\t', index=False)

    print('Sucessed to achieve gene annotation of "%s" from Ensembl biomart.' % _)


def search (name):
    msg = 'NotFound'
    name = name.replace(' ', '_').replace('-', '_')
 
    UA = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 ' + \
    '(KHTML, like Geko) chrom/61.0.3163.100'

    for i in ['www', 'metazoa', 'plants', 'fungi', 'bacteria', 'protists']:
        url = 'http://%s.ensembl.org/%s' % (i, name)
        
        query = requests.get(url,
        headers = {'User-Agent': UA, 'Referer': 'http://asia.ensembl.org'})

        # print('Search "%s" in http://%s.ensembl.org' % (name, i))
        if query.status_code != 200: continue
        msg = query.url; break

    print(msg)

    return(1 if msg == 'NotFound' else 0)


# arg1 = 'http://plants.ensembl.org/Glycine_max/Info/Index'
# arg1 = 'http://asia.ensembl.org/Mus_musculus/Info/Index'
cmd, arg1 = os.sys.argv[1:3]

if cmd == 'search':
    os.sys.exit( search (formatSpeciesName (arg1)))
elif cmd == 'getftp':
    netloc, path, loca = query(arg1)
    getftp (netloc, path, loca, arg1)
elif cmd == 'biomart':
    netloc, path, loca = query(arg1)
    biomart_anno(arg1, loca)
else:
    pass
