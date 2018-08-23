mkdir example

go run Pathway_download.go hsa

mv hsa00001.keg.gz example/

python3 Pathway_keg2tsv.py example/hsa00001.keg.gz example/hsa00001

Rscript pathway_classfication.r example/hsa00001.classfication.tsv example/hsa00001.classfication.pdf
