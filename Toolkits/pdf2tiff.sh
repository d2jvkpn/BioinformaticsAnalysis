#! /bin/bash

# author: d2jvkpn
# version: 0.3
# project: https://github.com/d2jvkpn/BioinformaticsAnalysis
# date: 2017-07-03
# description: convert pdf image to tiff format with default density 300.

if [ -z $(which convert) ]; then
   echo "Can't use command \"convert\", please install ImageMagick."
   exit
fi

usage=$(
cat << EOF
Convert pdf to tiff image, with default density 300, usage:
    sh  $(basename $0) [density_number]  <pdf_files>
EOF
)

if [ -z $1 ]; then echo "$usage"; exit; fi

test -f $1 && dens=300 || { dens=$1; shift; }

for i in $@; do
    if [[ ! "$i" == *".pdf" ]] || \
    [ ! -s $i ]; then echo "$i: Skip"; continue; fi

    convert -density $dens "$i" ${i%.pdf}.tiff 2>/dev/null

    test $? -eq 0 && sta="OK" || sta="Failed"
    echo "$i: $sta"
done
