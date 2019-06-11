#! /bin/bash

set -eu -o pipefail

# https://github.com/d2jvkpn/BioinformaticsAnalysis

# http://genomeview.org/manual/Bam2tdf
# http://genomeview.org/jenkins/bam2tdf-nightly/bam2tdf-14.zip
# java -jar bam2tdf.jar [-m <minimumMappingQuality>] <location of your bam file>

for i in $(ls *.sorted.bam); do
    java -jar bam2tdf-14/bam2tdf.jar $i &
done
