#! /bin/bash

python3 Species_txid.py

awk 'BEGIN{FS=OFS="\t"} {split($3,x,"[()]"); gsub(";", "; ", $4);
sub(" $", "", x[1]); print $1,$2,x[1],x[2],$4}' data/meta_organism.tsv |
awk 'BEGIN{FS=OFS="\t"; print "Entry", "code", "Scientific_name", "txid", \
"Common_name", " Taxonomic_classification"} NR==FNR{a[$1]=$2; next}
{$3=$3"\t"a[$3]; print}' data/Species2txid.tsv - > data/KEGG_organism.tsv
