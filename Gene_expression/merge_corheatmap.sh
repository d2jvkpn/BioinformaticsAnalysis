#! /bin/bash

set -eu
prefix=$1

# sudo add-apt-repository ppa:malteworld/ppa
# sudo apt update
# sudo apt install pdftk

for i in pearson spearman kendall; do
	test -s $prefix.corheatmap_$i.pdf || exit
done

pdftk $prefix.corheatmap_pearson.pdf $prefix.corheatmap_spearman.pdf \
$prefix.corheatmap_kendall.pdf cat output ${prefix}.corrheatmap.pdf

rm $prefix.corheatmap_pearson.pdf $prefix.corheatmap_spearman.pdf \
$prefix.corheatmap_kendall.pdf
