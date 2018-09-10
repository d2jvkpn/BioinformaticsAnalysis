# Update local data table for command "match":
PathwayTK  update

# Find match species in local data table:
PathwayTK  match  "Rhinopithecus roxellana"
PathwayTK  match  Rhinopithecus+roxellana

# Get organisms keg file:
PathwayTK  get  hsa mmu ath

# Convert keg format to tsv:
PathwayTK  totsv  hsa00001.keg.gz  hsa00001.keg.tsv
PathwayTK  totsv  hsa00001.keg  hsa00001.keg.tsv

# Get species keg and convert to tsv:
PathwayTK  species  Rhinopithecus+roxellana

sh gene2pathway.sh hsa00001.keg.tsv
