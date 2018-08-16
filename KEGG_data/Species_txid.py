import os, re
from bs4 import BeautifulSoup

html = 'data/catalog_Species.html'

with open(html, 'r') as f: HTML = f.read()

F = open('data/Species2txid.tsv', 'w')

for i in HTML.split('<tr align=center>'):
    ii = i.split('</tr>')[0]
    if '&category_type=species' not in ii: continue
    tds = BeautifulSoup(ii, 'html.parser').find_all('td')
    _ = tds[-2].select('a')
    txid = '' if _ == [] else tds[-2].select('a')[0].text
    F.write(tds[-3].select('a')[0].text + '\t' + txid + '\n')

f.close()
