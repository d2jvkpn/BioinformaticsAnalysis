#! /bin/bash
# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

set -eu -o pipefail

makeblastdb -in target.faa -input_type fasta -dbtype prot

blastx -query query.fa -db target.faa -evalue 1E-5 -strand plus \
-max_target_seqs 5 -num_threads 10 -outfmt 6 |
awk 'BEGIN{FS=OFS="\t"; print "qseqid", "sseqid", "pident", "length", \
"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"}
{print}' > blastx_target.tsv
