#! /bin/bash

set -eu -o pipefail

tsv=$1

awk 'BEGIN{FS=OFS="\t"} $2~"^PATH:"{sub("PATH:", "", $2); a[$1]=$2; next}
a[$1] && $2=="-"{print $4,a[$1]":"$10":"$12}' $tsv |
sort -u | awk 'BEGIN{FS=OFS="\t"; print "gene", "pathway"}
{print}' > gene2pathway.tsv

awk 'BEGIN{FS=OFS="\t"; print "pathway", "pathway_name", "A", "B"}
$2~"^PATH:"{sub("PATH:", "", $2);
print $2,$3,$7,$9; next}' $tsv > pathway.infor.tsv
