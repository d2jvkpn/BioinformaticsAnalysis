python3 Gffinfor_v0.6.py GCF_000001405.38_GRCh38.p12_genomic.gff.gz
#    NAME                                              COUNT

#    Seqences                                          594
#    source BestRefSeq                                 1375428
#    source BestRefSeq%2CGnomon                        11832
#    source Curated Genomic                            41242
#    source Gnomon                                     2221968
#    source RefSeq                                     42350
#    source tRNAscan-SE                                1764
#    type antisense_RNA                                22
#    type biological_region                            2233
#    type CAAT_signal                                  5
#    type CAGE_cluster                                 94
#    type cDNA_match                                   15196
#    ......


python3 Gffinfor_v0.6.py GCF_000001405.38_GRCh38.p12_genomic.gff.gz gene
#    TYPE  ATTRIBUTION   TOTAL  UNIQUE

#    gene  Dbxref        42955  38489
#    gene  description   31057  27411
#    gene  end_range     572    568
#    gene  exception     16     1
#    gene  gbkey         42955  1
#    gene  gene          42955  38474
#    gene  gene_biotype  42955  21
#    gene  gene_synonym  22683  19654
#    gene  ID            42955  42955
#    gene  Name          42955  38474
#    gene  partial       851    1
#    gene  start_range   487    463


python3 Gffinfor_v0.6.py GCF_000001405.38_GRCh38.p12_genomic.gff.gz \
gene,pseudogene gene,gene_biotype,description GeneID,HGNC |
awk '++x[$1]==1' > gene.infor.tsv


####
python3 NCBI_gff2gtf_v1.7.py GCF_000001405.38_GRCh38.p12_genomic.gff.gz raw_gene

python3 NCBI_gtf4stringtie.py raw_gene.gtf.gz > gene.gtf

awk 'BEGIN{FS=OFS="\t"} $9~"transcript_biotype \"mRNA\""{print}' gene.gtf > mRNA.gtf
