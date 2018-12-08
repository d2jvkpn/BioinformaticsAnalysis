#! /bin/bash

set -eu -o pipefail

USAGE='''
Convert KAAS output keg format to tsv, usage:
  $ sh KAAS_keg2tsv.sh  <input.keg>  [output.tsv]

author: d2jvkpn
version: 0.2
release: 2018-12-07
project: https://github.com/d2jvkpn/BioinformaticsAnalysis
lisense: GPLv3  (https://www.gnu.org/licenses/gpl-3.0.en.html)
'''

if [ $# -lt 1 ] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    echo "$USAGE"; exit 2
fi

keg=$1

if [ $# -eq 1 ]; then
    tsv=/dev/stdout
else
    tsv=$2
fi

eval cat $keg |
awk 'BEGIN{OFS="\t"; print "target_id", "KO_id", "KO_information", "C_entry", 
"C_id", "C_name", "B_id", "B_name", "A_id", "A_name", "BR1", "BR2"}
$1~"^A"{sub(" ", "\t", $0); split($0,x,"\t"); Ai=x[1]; An=x[2]; next}
$1=="B" && NR>2{sub("  ", "", $0); sub(" ", "\t", $0); split($0,x,"\t");
  Bi=x[1]; Bn=x[2]}
$1=="C"{k=0; if($NF!~"BR:" && $NF!~"PATH:") next; k=1; Ci="C"$2;
  sub("C    ..... ", "", $0); Pe=$NF; sub("[[]", "", Pe); sub("[]]", "", Pe);
  sub(" [[].*", "", $0); Cn=$0}
$1=="D" && k==1{BR1=""; BR2=""; split($0,fd,"\t");
  if(length(fd)>1) BR1=fs[2]; if(length(fd)==3) BR2=fd[3];
  sub("D      ", "", fd[1]); split(fd[1], x, "  "); split(x[1], y, "; "); 
  print y[1], y[2], x[2], Pe, Ci, Cn, Bi, Bn, Ai,An, BR1, BR2}' > $tsv
