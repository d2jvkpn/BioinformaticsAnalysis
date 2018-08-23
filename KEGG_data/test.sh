
mkdir example

code=hsa

go run Pathway_download.go $code

mv ${code}00001.keg.gz example/

python3 Pathway_keg2tsv.py example/${code}00001.keg.gz example/${code}00001

Rscript pathway_classfication.r example/${code}00001.classfication.tsv \
example/${code}00001.classfication.pdf
