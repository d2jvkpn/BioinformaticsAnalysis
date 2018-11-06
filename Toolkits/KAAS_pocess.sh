#! /bin/bash

set -eu -o pipefail

keg=$1
class=$2
classtsv=$(dirname $(which Pathway))/KEGG_data/LineageL2_PATH.tsv

test -s $keg || exit
test -s $classtsv || exit

Pathway tsv $keg |
awk -v c=$class 'BEGIN{FS=OFS="\t"} NR==FNR{split($1,x,"");
if(NR>1 && x[1]==c) a[$2]=1; next} /^PATH:/{sub("PATH:", "", $1);
if(a[$1]) print $2,$5}' $classtsv - | sort -u  >  gene2ko.txt
