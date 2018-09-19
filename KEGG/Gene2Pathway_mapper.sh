#! /bin/bash

set -eu -o pipefail

tsv=$1

cat $tsv | awk 'BEGIN{FS=OFS="\t"}
$1~"^#" && $2~"^PATH:"{sub("#", "", $1); a[$1]=$2;
  sub(".*:", "", $4); sub(".*:", "", $5);
  b[$1]=$3"\t"$4"\t"$5; next}
a[$1] {sub("PATH:", "", a[$1]); split($3,x,"; ");
if(length(x)==2) {g=x[1]} else {g=""};
sub("[[]EC:", "", $4); sub("[]]", "", $4);
split($4, x, " "); print $2, a[$1], g, x[1], $5, b[$1]}'
