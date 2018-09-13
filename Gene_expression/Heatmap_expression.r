library (pheatmap)

Args <- commandArgs(T)
if (length(Args) == 0 || Args[1] == "-h" || Args[1] == "--help") {
	print("Arguments:  <expression.tsv>  <output.pdf>")
	print("author: d2jvkpn")
	print("version: 0.1")
	print("release: 2018-09-13")
	print("project: https://github.com/d2jvkpn/BioinformaticsAnalysis")
	print("lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html")
	q()
}


hd <- read.delim (Args[1], check.names=FALSE, row.names=1)
PDF <- Args[2]

hd <- hd [apply(hd,1,var) > 0,]

palette <- colorRampPalette (c ("blue", "white", "red"))(n=299)
print (PDF)

pheatmap (log (hd + 1, 10),
    scale = ifelse (ncol(hd) > 2, "row", "column"),
    breaks=NA, color=palette, cluster_rows = TRUE,
    cluster_cols = ncol(hd) > 2, show_rownames = nrow(hd) <= 50,
    border_color=NA, filename=PDF,  heigt=10, width=6)
