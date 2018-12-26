#! /bin/bash
# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

set -eu -o pipefail

makeblastdb -in target.fa -input_type fasta -dbtype nucl

blastn -query query.fa -db target.fa \
-perc_identity 80 -gapopen 0 -evalue 1E-5 -strand plus \
-max_target_seqs 5 -num_threads 10 -outfmt 6 |
awk 'BEGIN{FS=OFS="\t"; print "qseqid", "sseqid", "pident", "length", \
"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"}
{print}' > blastn_target.tsv

# qseqid    Query Seq-id
# sseqid    Subject Seq-id
# pident    Percentage of identical matches
# length    Alignment length
# mismatch  Number of mismatches
# gapopen   Number of gap openings
# qstart    Start of alignment in query
# qend      End of alignment in query
# sstart    Start of alignment in subject
# send      End of alignment in subject
# evalue    Expect value
# bitscore  Bit score
