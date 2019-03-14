#! /bin/bash
# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

set -eu -o pipefail

makeblastdb -in target.fa -input_type fasta -dbtype nucl

evalue=1E-10
identity=80
out=blastn_target.tsv

blastn -query query.fa -db target.fa -perc_identity $identity -gapopen 0 \
-evalue $evalue -strand plus -max_target_seqs 5 -num_threads 10 -outfmt 6 |
awk 'BEGIN{FS=OFS="\t"; print "qseqid", "sseqid", "pident", "length", \
"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", \
"bitscore"} {print}' > $out

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
