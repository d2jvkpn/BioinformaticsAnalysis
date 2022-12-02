#! /usr/bin/env bash
set -eu -o pipefail
_wd=$(pwd)
_path=$(dirname $0 | xargs -i readlink -f {})

mkdir -p blastn

for subj in subject/*.fasta; do
    x=$(basename $subj | sed 's/.fasta$//')
    mkdir -p blastn/$x

    for fa in clean/*.fasta; do
        s=$(basename $fa | sed 's/.fasta$//')
        echo ">>>" $s
 
        blastn -query $fa -db $subj -word_size 7 -gapopen 0 -max_target_seqs 5 \
          -num_threads 2 -outfmt 5 |
          pigz -c > blastn/$x/blastn_$s.raw.xml.gz

        python3 scripts/blastn_part_xml2tsv.py blastn/$x/blastn_$s.raw.xml.gz \
          blastn/$x/blastn_$s.raw.tsv
    done
done
