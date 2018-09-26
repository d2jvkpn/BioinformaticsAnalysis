#! /bin/bash

set -eu -o pipefail

wd=$PWD
which Pathway
dp=$(dirname $(which Pathway))

cd $dp

awk -F "\t" 'NR>1{print $2}' KEGG_organism.tsv |
xargs -i -n 100 Pathway Get &> get_pathway.log

awk '/Failed/{print $NF}' get_pathway.log > get_pathway.failed

for i in $(ls *.keg.gz); do
    gunzip -c $i | grep -w -o "PATH:.*" |
    awk -v i=$(basename $i | sed 's/[0-9]*.keg.gz//') 'BEGIN{FS=OFS="\t"}
    {sub("PATH:", "", $1); sub("[]]", "", $1); sub(i, "ko", $1); print i, $1}'
done |
awk 'BEGIN{FS=OFS="\t"} NR==FNR{if(NR>1) a[$2]=$4; next}
a[$1]{if(++x[a[$1]"\t"$2]) print a[$1], $2}' KEGG_organism.tsv - |
sort | sed '1i Lineage\tPATH' | gzip -c > Lineage_PATH.tsv.gz

gunzip -c Lineage_PATH.tsv.gz | awk 'BEGIN{FS=OFS="\t"}
NR==1{print} NR>1{split($1,x,";");
if (++y[x[2]"\t"$2]==1) print x[2],$2}' > LineageL2_PATH.tsv

tar -cf Pathway_keg.tar *.keg.gz
rm *.keg.gz

# tar -xf Pathway_keg.tar hsa00001.keg.gz

cd $wd
