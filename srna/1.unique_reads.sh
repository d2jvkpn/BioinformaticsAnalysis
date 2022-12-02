#! /usr/bin/env bash
set -eu -o pipefail
_wd=$(pwd)
_path=$(dirname $0 | xargs -i readlink -f {})

mkdir -p clean
NC=4; NG=4

for fq in $(ls fastq/*.fq.gz); do
    s=$(basename $fq); s=${s%.fq.gz}

    pigz -dc $fq | awk 'NR%4==2{print $1}' |
      sort --parallel=$NC --buffer-size=${NG}G |
      uniq -c |
      awk 'BEGIN{n=0} {n++; printf ">t%07d %d\n%s\n", n, $1, $2}' > clean/$s.fasta

    awk 'BEGIN{OFS="\t"; print "id", "sequence", "length", "copy"}
      NR%2==1{sub(">", ""); i=$1; c=$2}
      NR%2==0{print i, $1, length($1), c}' clean/$s.fasta |
      pigz -c > clean/$s.fasta.tsv.gz

    awk 'BEGIN{FS=OFS="\t"} NR>1{a[$3]+=$4} END{print "length", "count";
      for (i in a) print i,a[i] | "sort -n"}' clean/$s.tsv > clean/$s.length.tsv
done

Rscript reads_length_dist.r clean/*.length.tsv
