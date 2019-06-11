from Bio import SeqIO

# https://github.com/d2jvkpn/BioinformaticsAnalysis

out = open("new.fa", "w")

for r in SeqIO.parse("input.fa", "fasta"):
    out.write(">" + r.id + "\n")
    out.write(str(r.seq)  + "\n")

out.close()
