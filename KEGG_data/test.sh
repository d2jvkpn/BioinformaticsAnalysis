
mkdir -p example

code=hsa

go run Pathway_download.go $code

mv ${code}00001.keg.gz example/

python3 Pathway_keg2tsv.py example/${code}00001.keg.gz example/${code}00001

Rscript pathway_classfication.r example/${code}00001.classfication.tsv \
example/${code}00001.classfication.pdf

####
python3 Pathway_find.py "Lokiarchaeum sp. GC14_75"  

# Entry: T04015
# code: loki
# Scientific_name: Lokiarchaeum sp. GC14_75
# txid: 
# Common_name: 
# Taxonomic_classification: Prokaryotes; Archaea; Lokiarchaeota; Lokiarchaeum
