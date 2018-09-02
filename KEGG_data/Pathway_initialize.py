#! /usr/bin/env python3

import os, re, subprocess
from bs4 import BeautifulSoup

cmd1 = '''
mkdir -p $SP/data

{
  wget -O meta_organism.tsv.part http://rest.kegg.jp/list/organism -o /dev/null &&
  mv meta_organism.tsv.part $SP/data/meta_organism.tsv || rm meta_organism.tsv.part
} &
p1=$!

{
  wget -O meta_pathway.tsv.part http://rest.kegg.jp/list/pathway -o /dev/null &&
  mv meta_pathway.tsv.part $SP/data/meta_pathway.tsv || rm meta_pathway.tsv.part
} &
p2=$!

{
  wget -O catalog_Genomes.html.part https://www.kegg.jp/kegg/catalog/org_list.html -o /dev/null &&
  mv catalog_Genomes.html.part $SP/data/catalog_Genomes.html || rm catalog_Genomes.html.part
} &
p3=$!

{
  wget -O catalog_Species.html.part https://www.kegg.jp/kegg/catalog/org_list4.html -o /dev/null &&
  mv catalog_Species.html.part $SP/data/catalog_Species.html || rm catalog_Species.html.part
} &
p4=$!

{
  wget -O catalog_Genus.html.part https://www.kegg.jp/kegg/catalog/org_list5.html -o /dev/null &&
  mv catalog_Genus.html.part $SP/data/catalog_Genus.html || rm catalog_Genus.html.part
}
p5=$!

{
  wget -O catalog_Viruses.html.part https://www.kegg.jp/kegg/catalog/org_list2.html -o /dev/null &&
  mv catalog_Viruses.html.part $SP/data/catalog_Viruses.html || rm catalog_Viruses.html.part
}
p6=$!

{
  wget -O catalog_Meta.html.part https://www.kegg.jp/kegg/catalog/org_list3.html -o /dev/null &&
  mv catalog_Meta.html.part $SP/data/catalog_Meta.html || rm catalog_Meta.html.part
} &
p7=$!

wait $p1 $p2 $p3 $p4 $p5 $p6 $p7
'''

cmd2 = '''
awk 'BEGIN{FS=OFS="\t"} {split($3,x,"[()]"); gsub(";", "; ", $4);
sub(" $", "", x[1]); print $1,$2,x[1],x[2],$4}' $SP/data/meta_organism.tsv |
awk 'BEGIN{FS=OFS="\t"; print "Entry", "code", "Scientific_name", "txid", \
"Common_name", "Taxonomic_classification"} NR==FNR{a[$1]=$2; next}
{$3=$3"\t"a[$3]; print}' $SP/data/species2txid.tsv - > KEGG_organism.tsv &&
mv KEGG_organism.tsv $SP/data/
'''

SP = os.path.dirname (os.path.realpath (__file__) )
status1, _ = subprocess.getstatusoutput (('SP=%s\n' % SP) + cmd1)
if status1 != 0: os.sys.exit ("Error, failed to download html.")

with open (SP + '/data/catalog_Species.html', 'r') as f: HTML = f.read ()
F = open (SP + '/data/species2txid.tsv', 'w')

for i in HTML.split ('<tr align=center>'):
    ii = i.split ('</tr>')[0]
    if '&category_type=species' not in ii: continue
    tds = BeautifulSoup (ii, 'html.parser').find_all('td')
    _ = tds[-2].select ('a')
    txid = '' if _ == [] else tds[-2].select('a')[0].text
    F.write (tds[-3].select ('a')[0].text + '\t' + txid + '\n')

F.close ()

status2, _ = subprocess.getstatusoutput (('SP=%s\n' % SP) + cmd2)

if status2 != 0: os.sys.exit ("Error: failed to parse data")
