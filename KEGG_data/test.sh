# Update local data table for command "match":
Pathway  update

# Find match species in local data table:
Pathway  match  "Rhinopithecus roxellana"
Pathway  match  Rhinopithecus+roxellana

# Get organisms keg file:
Pathway  get  hsa mmu ath

# Convert keg format to tsv:
Pathway  totsv  hsa00001.keg.gz  hsa00001.keg.tsv
Pathway  totsv  hsa00001.keg  hsa00001.keg.tsv

# Get species keg and convert to tsv:
Pathway  species  Rhinopithecus+roxellana
