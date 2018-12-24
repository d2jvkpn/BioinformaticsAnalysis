#! /bin/bash
set -eu -o pipefail


Rscript Enrichment_phyper_v0.2.r data/target_gene.tsv data/gene2GO.tsv \
2 data/GO_enrichment_target_gene.tsv

Rscript Enrichment_phyper_v0.2.r data/target_gene.tsv data/gene2pathway.tsv \
3 data/Pathway_enrichment_target_gene.tsv
