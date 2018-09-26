#! /bin/bash

# author: d2jvkpn
# date: 2017-07-28
# version: 1.08
# description: convert bam to cram, sam to cram, and cram to bam, make sure you
# have samtools 1.0 or higher.

usage=$(
cat <<EOF
Convert bam to cram, sam to cram, and cram to bam, USAGE:
    sh  $0  <bam/sam/cram files>  -f  <genome_fasta>  -p [threads]

    Note: Make sure you have samtools 1.0 or higher.
EOF
)

if [ -z $1 ] || [[ $1 == "-h" ]] || [[ $1 == "--help" ]]; then
    echo "$usage"; exit 0
fi

test -z $(which samtools) && { echo "Samtools is not available, quit."; exit; }

inputs=""
for i in $@; do
    test -s $i && { inputs="$inputs $i"; shift; } || break
    if [[ "$i" == "-f" ]] || [[ "$i" == "-p" ]]; then break; fi
done

while getopts "f:p:" arg; do
    case $arg in
        f) genome_fa=$OPTARG;;
        p) threads="-@ $OPTARG";;
        *) echo "$usage"; exit 0;;
    esac
done

test -s $genome_fa || { echo "File \"$genome_fa\" not exists."; exit; }

{
echo "StartTime: $(date +"%Y-%m-%d %H:%M:%S %z")"

for i in $inputs; do
    test -s $i || { echo "$i: NotExists"; continue; }

    if [[ $i == *".bam" ]] || [ $i == *".sam" ]; then
        o=${i%.bam}; o=${o%.sam}.cram; bC="-C"
    elif [[ $i == *".cram" ]]; then
        o=${i%.cram}.bam; bC="-b"
    else
        echo "$i: Error"; continue
    fi

    test -s $o && { echo "$i: Skip"; continue; }

    >&2 echo "Converting $i..."
    samtools view $threads -T $genome_fa $i $bC -o ${o}...

    test $? -eq 0 && { mv ${o}... $o; echo "$i: OK"; } ||
    { rm $o; echo "$i: Failed"; }
done

echo "EndTime: $(date +"%Y-%m-%d %H:%M:%S %z")"
} > bam2cram_$$.log
