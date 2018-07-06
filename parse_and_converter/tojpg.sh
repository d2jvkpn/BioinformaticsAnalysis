#! /bin/bash

# project: https://github.com/d2jvkpn/BioinformaticsAnalysis

if [ -z $(which convert) ]; then
   echo "Can't use command \"convert\", please install imagemagick."
   exit
fi

usage=$(
cat << EOF
Convert picture(s) to jpg format, with density 150 and quality 100.
Usage:
    sh  $(basename $0)  <input_files>
EOF
)

if [ $# -eq 0 ]; then echo "$usage"; exit; fi

for i in $@; do
    if [ ! -f $i ] || [ ! -s $i ]; then
        echo "$i: Not Exists."; continue
    fi
    
    convert -verbose -density 150 -trim  $i -quality 100 \
    -flatten -sharpen 0x1.0 ${i%.*}.jpg 2>/dev/null
    
    test $? -eq 0 && sta="OK" || sta="Failed"
    echo $i": "$sta

done
