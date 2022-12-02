#! /usr/bin/env bash
set -eu -o pipefail
_wd=$(pwd)
_path=$(dirname $0 | xargs -i readlink -f {})

fasta=$1

awk '/^>/{id=$1; sub(">", "", id)} {print $0 > id".fasta"}' $fasta
