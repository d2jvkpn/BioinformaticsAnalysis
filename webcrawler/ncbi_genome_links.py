#! /usr/bin/env python3

__author__ = 'd2jvkpn'
__version__ = '0.1'
__release__ = '2018-06-21'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

# export LANG='en_US.UTF-8'
# export LC_ALL='en_US.UTF-8'

import requests, os, time
from bs4 import BeautifulSoup

if len(os.sys.argv) == 1 or os.sys.argv[1] in ['-h', '--help']:
    print(''' 
Usage: python3 NCBI_genome_links.py "species scientific name" | "genome link"
e.g. Homo sapiens.''')
    os.sys.exit(0)


url = ' '.join(os.sys.argv[1:])

if url.count('ncbi.nlm.nih.gov') == 0:
    url = 'https://www.ncbi.nlm.nih.gov/genome/?term=' + url.replace(' ', '+')

print('# URL: ' + url, end='\n\n')

r = requests.get(url)

html_soup = BeautifulSoup(r.text, 'html.parser')

print('# Time: ' + time.strftime ('%Y-%m-%d %H:%M:%S %z'), end='\n\n')

Lineage = html_soup.find('span', class_='GenomeLineage').find_all('a')
print('# Lineage name: ' + ', '.join([i.text  for i in Lineage[0::2]]), end='\n\n')

_ = [i.get('href').split('/')[-1] for i in Lineage[0::2]]
print('# Lineage txid: ' + ', '.join(_), end='\n\n')

RGS = html_soup.find('div', class_='refgenome_sensor')

if len(RGS.select('b')) == 3:
    print('# Reference genome:', RGS.select('a')[0].text, end = '\n\n')

_ = [i.find_all('a') for i in RGS.find_all ('span', class_='shifted')]

# response = requests.get('http://www.example.com/image.jpg', stream=True)

for i in sum(_, []):
    href = i.get('href')
    if href.startswith('ftp://') and os.path.basename(href) != '':
        print('# download ' + i.text, href, sep='\n', end='\n\n')
