# GenomicProcess
tools for genomic data analysis (process)

**[ncbi_gff3_to_ensembl_gtf.py]**

  Convert ncbi gff3 format to ensembl gtf, usage:

    python3 ncbi_gff3_to_ensembl_gtf.py <input.gff3> <output.gtf>

  Note: transcription mapping file "genomic.transcription.tsv" will be created.

<br></br>

**[gff3infor]**

Gff3 file statistic and information extraction, usage:

    sh  gff3infor  <input.gff3>  [feature-type]  [attributions]

Note: gawk4.0 or higher is required to support urldecode
