#! /bin/bash

mkdir -p data

{
  wget -O meta_organism.tsv.part http://rest.kegg.jp/list/organism &&
  mv meta_organism.tsv.part data/meta_organism.tsv || rm meta_organism.tsv.part
} &

{
  wget -O meta_pathway.tsv.part http://rest.kegg.jp/list/pathway &&
  mv meta_pathway.tsv.part data/meta_pathway.tsv || rm meta_pathway.tsv.part
} &


{
  wget -O catalog_Genomes.html.part https://www.kegg.jp/kegg/catalog/org_list.html &&
  mv catalog_Genomes.html.part data/catalog_Genomes.html || rm catalog_Genomes.html.part
} &

{
  wget -O catalog_Species.html.part https://www.kegg.jp/kegg/catalog/org_list4.html &&
  mv catalog_Species.html.part data/catalog_Species.html || rm catalog_Species.html.part
} &

{
  wget -O catalog_Genus.html.part https://www.kegg.jp/kegg/catalog/org_list5.html &&
  mv catalog_Genus.html.part data/catalog_Genus.html || rm catalog_Genus.html.part
}

{
  wget -O catalog_Viruses.html.part https://www.kegg.jp/kegg/catalog/org_list2.html &&
  mv catalog_Viruses.html.part data/catalog_Viruses.html || rm catalog_Viruses.html.part
}

{
  wget -O catalog_Meta.html.part https://www.kegg.jp/kegg/catalog/org_list3.html &&
  mv catalog_Meta.html.part data/catalog_Meta.html || rm catalog_Meta.html.part
} &

wait

python3 Species_txid.py

awk 'BEGIN{FS=OFS="\t"} {split($3,x,"[()]"); gsub(";", "; ", $4);
sub(" $", "", x[1]); print $1,$2,x[1],x[2],$4}' data/meta_organism.tsv |
awk 'BEGIN{FS=OFS="\t"; print "Entry", "code", "Scientific_name", "txid", \
"Common_name", " Taxonomic_classification"} NR==FNR{a[$1]=$2; next}
{$3=$3"\t"a[$3]; print}' data/Species2txid.tsv - > data/KEGG_organism.tsv
