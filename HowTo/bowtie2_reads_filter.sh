#! /bin/bash
set -eu -o pipefail

bowtie2-build -f target.fa target

mkdir -p log

bowtie2 -k30 -t -p 15 -x target -1 $i.R1.fastq.gz -2 $i.R2.fastq.gz \
--un-conc ./un-conc.$i --al-conc ./al-conc.$i -S /dev/null 2> log/bowtie2.$i.log

mv ./un-conc.$i.1 un-conc_$i.R1.fastq
mv ./un-conc.$i.2 un-conc_$i.R1.fastq
mv ./al-conc-mate.$i.1 al-conc_$i.R1.fastq
mv ./al-conc-mate.$i.2 al-conc_$i.R2.fastq

pigz un-conc_$i.R1.fastq un-conc_$i.R1.fastq \
al-conc_$i.R1.fastq al-conc_$i.R2.fastq
