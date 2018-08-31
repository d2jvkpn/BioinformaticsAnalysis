python3 NCBI_gff2gtf_v1.7.py GCF_000001405.38_GRCh38.p12_genomic.gff.gz raw_gene

python3 NCBI_gtf4stringtie.py raw_gene.gtf.gz > gene.gtf

awk 'BEGIN{FS=OFS="\t"} $9~"transcript_biotype \"mRNA\""{print}' gene.gtf > mRNA.gtf
