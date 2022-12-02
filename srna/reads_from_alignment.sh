#! /usr/bin/env bash
set -eu -o pipefail
_wd=$(pwd)
_path=$(dirname $0 | xargs -i readlink -f {})

mkdir clean

for f in $(ls blastn_ref/*/*.align.tsv); do
    s=$(basename $f | sed 's/^blastn_//; s/.align.tsv$//')
    echo ">>>" $s

    awk 'BEGIN{FS="\t"} NR>1{print ">"$1" "$13; print $2}' $f > clean/$s.fasta

    awk 'BEGIN{FS=OFS="\t"; print "id", "sequence", "length", "copy"}
      NR>1{print $1, $2, $8, $13}' $f |
      pigz -c > clean/$s.fasta.tsv.gz
done
