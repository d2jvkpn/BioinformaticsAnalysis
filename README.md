## parse_and_converter

Tools for Genomic Data Analysis

**[ncbi_gff3_to_ensembl_gtf.py]**

  Convert ncbi gff3 (gzipped) format to ensembl gzipped gtf (be compatible with softwares like stringtie), usage:

    python3  ncbi_gff3_to_ensembl_gtf.py  <input.gff3.gz>  <outputPrefix>


**[gffinfor.py]**

Extract attribution(s) of feature(s) from gtf or gff3 (gzipped file is supported), usage:

    python3  gtfinfor.py  <input.gtf/input.gff3>  <feature1,feature2...>  <attribution1,attribution2...>

Output: stdout, tsv format.



**[gff3infor.sh]**

Gff3 file statistic and information extraction, usage:

    sh  gff3infor.sh  <input.gff3>  [feature-type]  [attribution1,attribution2...]

Note: gawk4.0 or higher is required to support urldecode

***
## webcrawler

**[ncbi_genome_links.py]**
extract organism information and download links (fna, faa, gff, gbff) from genome webpage.

**[SRR_ncbi_bioproject.py]**
extact SRX and SRR data from ncbi bioproject webpage.

**[SRR_download]**
download sequecing fastq file by read tsv creatd by SRR_ncbi_bioproject.py.
