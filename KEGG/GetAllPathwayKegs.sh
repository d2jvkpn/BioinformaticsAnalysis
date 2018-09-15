#! /bin/bash

set -eu -o pipefail

wd=$PWD
which Pathway
dp=$(dirname $(which Pathway))

cd $dp

awk -F "\t" 'NR>1{print $2}' KEGG_organism.tsv |
xargs -i -n 100 Pathway Get &> get_pathway.log

awk '/Failed/{print $NF}' get_pathway.log > get_pathway.failed

mkdir Pathway_kegs && mv *.keg.gz Pathway_kegs/

cd $wd
