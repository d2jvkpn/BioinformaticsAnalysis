#! /bin/bash
# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

set -eu -o pipefail

bowtie2-build -f target.fa target

mkdir -p log

bowtie2 -k30 -t -p 15 -x target -1 $i.R1.fastq.gz -2 $i.R2.fastq.gz \
--un-conc ./un-conc.$i --al-conc ./al-conc.$i -S /dev/null 2> log/bowtie2.$i.log

mv ./un-conc.1.$i un-conc_$i.R1.fastq
mv ./un-conc.2.$i un-conc_$i.R2.fastq
mv ./al-conc-mate.1.$i al-conc_$i.R1.fastq
mv ./al-conc-mate.2.$i al-conc_$i.R2.fastq

pigz un-conc_$i.R1.fastq un-conc_$i.R2.fastq \
al-conc_$i.R1.fastq al-conc_$i.R2.fastq
