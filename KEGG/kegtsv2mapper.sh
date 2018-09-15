#! /bin/bash

set -eu -o pipefail

tsv=$1

awk 'BEGIN{FS=OFS="\t"} $2~"^PATH:"{a[$1]=$2; b[$1]=$3"\t"$7"\t"$9; next}
a[$1] {sub("PATH:", "", a[$1]); split($5,x,"; ");
if(length(x)==2) {g=x[1]} else {g=""};
sub("[[]EC:", "", $11); sub("[]]", "", $11);
split($1, x, " "); print $4, a[$1], g, x[1], $11, b[$1]}' $tsv
