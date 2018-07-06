#! /bin/env python3

__author__ = 'd2jvkpn'
__version__ = '0.1'
__release__ = '2018-06-21'
__project__ = 'https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__ = 'GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'


import requests, os
from bs4 import BeautifulSoup

url = os.sys.argv[1]

if url.startswith('https://'):
    ur = 'https://www.ncbi.nlm.nih.gov/bioproject/' + url

r = requests.get(url)

html_soup = BeautifulSoup(r.text, 'html.parser')

sra = ''

for i in html_soup.find_all('td', class_='alignRight'):
    _ = i.select('a')[0].get('href')
    if _.startswith('/sra?linkname='): sra = _

if sra == '':
    os.sys.exit('Error: failed to find SRA link.')

r = requests.get( 'https://www.ncbi.nlm.nih.gov' + sra)

html_soup = BeautifulSoup(r.text, 'html.parser')

srx = {}

for i in html_soup.find_all('p', class_ = 'title'):
    srx[i.text] = 'https://www.ncbi.nlm.nih.gov' + i.select('a')[0].get('href')

_ = 'SRR_of_BioProject_' + os.path.basename(url.strip('/')) + '.tsv'

with open(_, 'w') as f:
    f.write('\t'.join(['SRA', 'SRX', 'SRR', 'Layout']) + '\n')

    for i in srx:
        r1 = requests.get(srx[i])
        html_soup = BeautifulSoup(r1.text, 'html.parser')
        SRR = html_soup.find('table').find('a').text
        _ = html_soup.find_all('div', class_='expand-body')[-1]
        Layout = _.find_all('div')[-1].select('span')[0].text

        f.write('\t'.join([i, srx[i], SRR, Layout]) + '\n')
