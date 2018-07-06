#! /bin/bash

# project: https://github.com/d2jvkpn/BioinformaticsAnalysis

i=$1

convert -density 300 -depth 8 -quality 85 $i ${i%.*}.png 2> /dev/null

test $? -eq 0 && echo "OK: $i" || echo "Failed: $i"
