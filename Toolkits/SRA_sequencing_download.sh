#! /bin/bash

__author__='d2jvkpn'
__version__='0.1'
__release__='2018-06-21'
__project__='https://github.com/d2jvkpn/BioinformaticsAnalysis'
__license__='GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html)'

mkdir -p log

SSR=$1

grep -v "^#" $SSR | awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$3,$4}' |
while read a b c; do
  mkdir -p $a
  if [[ "$c" == "PAIRED" ]]; then
    fastq-dump ---gzip --split-3 $b --outdir $a &> log/download_$a.log &
  else
    fastq-dump ---gzip $b --outdir $a &> log/download_$a.log &
  fi
done

#cache-mgr --report
#cache-mgr --clear
