Rscript gene_expression_group.r example/gene_TPM.tsv example/group/gene_TPM 0.01 \
TPM example/group.tsv

sh merge_corheatmap.sh example/group/gene_TPM


Rscript gene_expression_singles.r example/gene_TPM.tsv example/singles/gene_TPM 0.01 TPM

sh merge_corheatmap.sh example/singles/gene_TPM
