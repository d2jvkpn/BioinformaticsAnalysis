library (pheatmap)

Args <- commandArgs (T)
if (length (Args) < 2 || length (Args) > 4 ) {
	print ("Arguments:  <expression.tsv>  <output.pdf>  [pdf_heigh] [pdf_width]")
	print ("author: d2jvkpn")
	print ("version: 0.2")
	print ("release: 2018-10-28")
	print ("project: https://github.com/d2jvkpn/BioinformaticsAnalysis")
	print ("lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html")
	q ()
}


hd <- read.delim (Args[1], check.names=FALSE, row.names=1)
PDF <- Args[2]

if (length (Args) >= 3) { heigh = as.integer(Args[3]) } else { heigh = 10 }
if (length (Args) == 4) { width = as.integer(Args[4]) } else { width = 6 }

hd <- hd [apply(hd,1,var) > 0,]

palette <- colorRampPalette (c ("blue", "white", "red"))(n=299)
print (PDF)

pheatmap (log (hd + 1, 10),
    scale = ifelse (ncol(hd) > 2, "row", "column"),
    breaks = NA, color = palette, cluster_rows = TRUE,
    cluster_cols = ncol(hd) > 2, show_rownames = nrow(hd) <= 50,
    border_color = NA, filename = PDF, heigh = heigh, width = width)
