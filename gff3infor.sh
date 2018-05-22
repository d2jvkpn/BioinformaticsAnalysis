#! /bin/bash

# author: d2jvkpn
# date: 2017-08-11
# version: 1.6.6
# project: https://github.com/d2jvkpn/GenomicProcess
# description: extract information from gff3 file , gawk 4.0 or higher is
#+ required to support urldecode.

if [ -z $1 ]; then
    echo "Gff3 file statistic and information extraction, usage:"
    echo "  sh  $0  <input.gff3>  [feature-type]  [attributions]"
    echo "  Note: gawk4.0 or higher is required to support urldecode."
    exit 0
elif [[ ! $1 == "-" ]] && [ ! -s $1 ]; then
    echo "File not exists: \"$1\"."; exit
elif [[ $1 == *".gz" ]]; then
    read_gff3="gunzip -c $1"
else
    read_gff3="cat $1"
fi

if [[ $1 == *".gtf" || $1 == *".gtf.gz" ]]; then
   read_gff3="$read_gff3 | awk 'BEGIN{FS=OFS=\"\\t\"}
    { sub(\"\\\";\$\", \"\", \$9); gsub(\" \\\"\", \"=\", \$9);
    gsub(\"\\\"; \", \";\", \$9); print}' "
fi

case $# in
    1) ## GFF column 1, 2, 3 count
        eval $read_gff3 | gawk 'BEGIN{FS=OFS="\t"}
        $1!~"^#"{f1[$1]=1; f2[$2]++; f3[$3]++}
        END{ print "SequenceCount:", length(f1);
        for(i in f2) print "Source:", i, f2[i];
        for(i in f3) print "Feature:", i, f3[i] | "sort -k3,3nr" }' |
        column -t -s $'\t' | sed 's/^/    /';;

    2) ## Attributions (column 9) count of selected feature
        fea=$2

        echo "Attribution(s) of \"$fea\" (name, total, unique):"

        eval $read_gff3 |
        gawk -v f3=$fea -F "\t" '($3==f3){gsub(/;/,"\t", $9); print $9}' |
        gawk 'BEGIN{FS=OFS="\t"}{for(i=1; i<=NF; i++) {split($i, x, "=");
        k=x[1]; a[k]++; if(a[k] && !b[$i]) c[k]++; b[$i]=1}}
        END{for(i in a) print i, a[i], c[i] | "sort -k2,2nr -k3,3nr | \
        column -t"}' | awk '$0!~":"{$0="    "$0} {print}' ;;

    *) ## Extract selected attributions value of selected feature
        fea=$2; shift; shift; atrr=$@

        echo $atrr | sed 's/ /\t/g'

        eval $read_gff3 | 
        gawk -v f3=$fea -F "\t" 'BEGIN{split(f3, a, ",");
        for(i in a) b[a[i]]=i} (b[$3]){ print $9 }' |
        gawk -v ks="$atrr" 'BEGIN{FS=";"; OFS="\t"; split(ks, k, " ")}
        {for(i=1; i<=NF; i++){split($i,x,"="); a[x[1]]=x[2]}; \
        info=a[k[1]]; for(i=2; i<=length(k); i++) info=info"\t"a[k[i]];
        print info; split("", a, ":")}' |
        gawk -niord '{printf RT ? $0chr("0x"substr(RT,2)) : $0}' RS=%.. ;;
esac
