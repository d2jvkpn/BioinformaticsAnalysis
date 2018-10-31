#! /bin/bash

set -eu -o pipefail

gunzip -c PlantTFDB-all_TF_pep.fas.gz | grep ">" |
sed 's/>//; s/ /\t/; s/|/\t/g' |
awk 'BEGIN{FS=OFS="\t"; print "GeneID", "TF_protein", "TF_name", 
"TF_family"} $2=="Glycine max"{k=toupper($1); sub("[.]", "_", k);
sub("[.].*", "", k); print k,$1,$3,$4}' > GLYMA0.tf.tsv

python3 zipDF.py
# output is GLYMA0.tf.tsv
