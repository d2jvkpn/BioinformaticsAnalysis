#! /bin/bash

set -eu -o pipefail

export ShowRownames="TRUE"

Rscript Heatmap_expression.r target_gene.TPM.tsv target_gene.TPM.heatmap.pdf \
"Target gene heatmap"
