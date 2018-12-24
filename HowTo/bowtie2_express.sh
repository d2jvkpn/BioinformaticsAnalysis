#! /bin/bash
set -eu -o pipefail

bowtie2-build -f target.fa target

mkdir -p express_$i/ log/

bowtie2 -k30 -t -p 15 -x target -1 $i.R1.fastq.gz -2 $i.R2.fastq.gz \
2> log/bowtie2.$i.log |
express --no-update-check target.fa --rf-stranded -o express_$i/

cut -f2,8,11 express_$i/results.xprs > express_$i/$i.express.tsv
