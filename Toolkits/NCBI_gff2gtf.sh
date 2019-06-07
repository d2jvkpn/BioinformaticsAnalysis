#! /bin/bash

set -eu -o pipefail

gffread -EFG NCBI.gff3 -T -o genomic.gtf
# missing gene and transcript records

Gffinfor genomic.gtf exon \
transcript_id,gbkey,Name,product,gene_id,gene_biotype,gene_name |
awk 'BEGIN{FS=OFS="\t"} ++x[$1]==1{print}' > transcript.infor.tsv
