#! /bin/bash
# https://github.com/d2jvkpn/BioinformaticsAnalysis

set -eu -o pipefail

USAGE='''circos_io.sh usage:
  $ circos_io.sh  <circos.conf>  <output_prefix>
'''

if [ $# -ne 2 ]; then
	echo "$USAGE"
	exit 2
fi

test -z $(which circos) && { echo "no circos installed"; exit; }
test -f $1 || { echo "conf file not exists: $1"; exit; }

conf=$(readlink -f $1)
prefix=$2
tmp=$(mktemp -u -p ./)
echo "circos temporary directory: $tmp"

mkdir -p $tmp && cd $tmp
circos -conf $conf -silent
cd ../

if [[ "$prefix" == *"/" ]] || [ -d $prefix ]; then
	prefix=$prefix"/circos"
elif [[ "$prefix" != *"circos" ]]; then
	prefix=$prefix".circos"
fi

mkdir -p $(dirname $prefix)
mv $tmp/circos.png $prefix.png
mv $tmp/circos.svg $prefix.svg
rm -r $tmp
echo "saved $prefix.png $prefix.svg"
