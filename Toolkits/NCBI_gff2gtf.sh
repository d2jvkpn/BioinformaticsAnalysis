#! /bin/bash

# https://github.com/d2jvkpn/BioinformaticsAnalysis

set -eu -o pipefail

gffread -EFG NCBI.gff3 -T -o genomic.gtf
# missing gene and transcript records

Gffinfor genomic.gtf exon \
transcript_id,gbkey,Name,product,gene_id,gene_biotype,gene_name,description  \
gbkey:transcript_biotype,Name:transcript_name,description:gene_description |
awk 'BEGIN{FS=OFS="\t"} ++x[$1]==1{print}' > transcript.infor.tsv
