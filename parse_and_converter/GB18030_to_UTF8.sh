#! /bin/bash

__project__='https://github.com/d2jvkpn/BioinformaticsAnalysis'

for i in $@; do
  iconv -t UTF-8 -f GB18030 $i > ...tmp.txt
  if [ $? -eq 0 ]; then
    mv ...tmp.txt $i; echo "Done: $i"
  else
    rm ...tmp.txt; echo "Failed: $i"
  fi
done
