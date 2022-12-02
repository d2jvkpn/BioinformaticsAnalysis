#! /usr/bin/env bash
set -eu -o pipefail
_wd=$(pwd)
_path=$(dirname $0 | xargs -i readlink -f {})

# makeblastdb -in $subj -input_type fasta -dbtype nucl -out $subj

subj=$1
subj_name=$(basename $subj | sed 's/\.fasta$//; s/\.fa$//')

# sample name
for fa in $(ls clean/*.fasta); do
    s=$(basename $fa | sed 's/\.fasta$//; s/\.fa$//')

    echo ">>>" $s
    mkdir -p blastn/$subj_name

    blastn -query $fa -db $subj -perc_identity 100 -word_size 7 \
      -gapopen 0 -max_target_seqs 5 -num_threads 2 -outfmt 6 -evalue 0.1 |
      awk 'BEGIN{FS=OFS="\t"; print "qseqid", "sseqid", "pident", "length",
      "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue",
      "bitscore"} {print}' > blastn/$subj_name/blastn_$s.raw.tsv

    pigz blastn/$subj_name/blastn_$s.raw.tsv
done
