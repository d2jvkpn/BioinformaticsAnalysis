#! /usr/bin/env bash
set -eu -o pipefail

input=$1
perc=$2  # 50, 60, 70...

total=$(pigz -dc $input | awk '{n++} END{print n/4}')
number=$(($total*$perc/100))

echo "extract $perc($number/$total) reads from "

pigz -dc $input |
awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' |
awk -v k=$number 'BEGIN{srand(systime() + PROCINFO["pid"]);} {s=x++<k?x-1:int(rand()*x);
if(s<k) R[s]=$0} END{for(i in R) print R[i]}' | 
awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4}'
