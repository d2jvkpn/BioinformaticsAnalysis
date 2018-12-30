#! /bin/bash
# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

set -eu -o pipefail

# A program to identify plant transcription factors (TFs), transcriptional regulators (TRs) 
# and protein kinases (PKs) from protein or nucleotide sequences and then classify individual 
# TFs, TRs and PKs into different gene families. iTAK currently provides both online and 
# standalone versions. 

# https://github.com/kentnf/iTAK/archive/v1.7a.tar.gz
# ftp://itak.feilab.net/pub/program/itak/database/
# ftp://itak.feilab.net/pub/program/itak/old/iTAK-1.7.tar.gz

iTAK.pl -f 3F -p 20 target.fa &> iTAK_target.log &


# tf_sequence.fasta: sequences of all identified TFs/TRs (FASTA format).
# tf_classification.txt: classification of all identified TFs/TRs. A tab-delimited txt file containing sequence IDs and their families.
# tf_alignment.txt: A tab-delimited txt file containing alignments of all identified TFs/TRs to the protein domain database.
# pk_sequence.fasta: sequences of all identified PKs (FASTA format).
# Shiu_classification.txt: classification of all identified protein kinases. A tab-delimited txt file containing sequence IDs and their corresponding protein kinase families.
# Shiu_alignment.txt: A tab-delimited txt file containing alignments of all identified protein kinases to the protein domain database.
