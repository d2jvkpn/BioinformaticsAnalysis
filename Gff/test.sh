Gffinfor GCF_000001405.38_GRCh38.p12_genomic.gff.gz
#    NAME                                                 COUNT
#    Sequences                                            594
#    source: BestRefSeq                                   1375428
#    source: BestRefSeq%2CGnomon                          11832
#    source: Curated Genomic                              41242
#    source: Gnomon                                       2221968
#    source: RefSeq                                       42350
#    source: tRNAscan-SE                                  1764
#    type: antisense_RNA                                  22
#    type: biological_region                              2233
#    type: C_gene_segment                                 43
#    type: CAAT_signal                                    5
#    type: CAGE_cluster                                   94
#    type: cDNA_match                                     15196
#    type: CDS                                            1461835
#    ...


Gffinfor GCF_000001405.38_GRCh38.p12_genomic.gff.gz gene,mRNA
#    TYPE    ATTRIBUTION       TOTAL     UNIQUE
#    gene    Dbxref            42954     38488
#    gene    description       31056     27410
#    gene    end_range         571       567
#    gene    exception         15        1
#    gene    gbkey             42954     1
#    gene    gene              42954     38473
#    gene    gene_biotype      42954     21
#    gene    gene_synonym      22682     19653
#    gene    ID                42954     42954
#    gene    Name              42954     38473
#    gene    partial           850       1
#    gene    start_range       486       462
#    mRNA    Dbxref            119298    113606
#    mRNA    end_range         419       224
#    mRNA    exception         6877      2
#    mRNA    gbkey             119298    1
#    mRNA    gene              119298    20190
#    mRNA    ID                119298    119298
#    mRNA    inference         6869      4925
#    mRNA    model_evidence    63554     58998
#    mRNA    Name              119298    113606
#    mRNA    Note              6865      488
#    mRNA    Parent            119298    22932
#    mRNA    partial           694       1
#    mRNA    product           119298    113210
#    mRNA    start_range       283       147
#    mRNA    transcript_id     119298    113606

Gffinfor GCF_000001405.38_GRCh38.p12_genomic.gff.gz \
gene,pseudogene gene,gene_biotype,description GeneID,HGNC |
awk '++x[$1]==1' > gene.infor.tsv


python3 NCBI_gff2gtf_v1.7.py GCF_000001405.38_GRCh38.p12_genomic.gff.gz raw_gene
# output: 

python3 NCBI_gtf4stringtie.py raw_gene.gtf.gz > gene.gtf

awk 'BEGIN{FS=OFS="\t"} $9~"transcript_biotype \"mRNA\""{print}' gene.gtf > mRNA.gtf
