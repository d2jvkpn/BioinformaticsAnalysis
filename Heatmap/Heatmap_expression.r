library (pheatmap)

Args <- commandArgs (T)
if (length (Args) != 3 ) {
	print ("Arguments:  <expression.tsv>  <output.pdf> <title>")
    print ("  read environment variables: PDFheigh, PDFwidth, COLOR1, COLOR2,")
    print ("    COLOR3, PDFwidth, COLOR1, COLOR2, COLOR3, ShowRownames")
	print ("author: d2jvkpn")
	print ("version: 0.3")
	print ("release: 2018-12-05")
	print ("project: https://github.com/d2jvkpn/BioinformaticsAnalysis")
	print ("lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html")
	q ("no", 2, FALSE)
}


hd <- read.delim (Args[1], check.names=FALSE, row.names=1)
PDF <- Args[2]
title <- Args[3]

PDFheigh <- as.numeric(Sys.getenv(x="PDFheigh", unset="8"))
PDFwidth <- as.numeric(Sys.getenv(x="PDFwidth", unset="6"))
COLOR1 <- Sys.getenv(x="COLOR1", unset="blue")
COLOR2 <- Sys.getenv(x="COLOR2", unset="white")
COLOR3 <- Sys.getenv(x="COLOR3", unset="white")
ShowRownames <- Sys.getenv(x="ShowRownames", unset="")

if (nrow(hd) <= 50 & ShowRownames == "") {
	ShowRownames == "TRUE"
}

if (ShowRownames == "") {
    ShowRownames = "FALSE"
}

hd <- hd [apply(hd,1,var) > 0,]

palette <- colorRampPalette (c (COLOR1, COLOR2, COLOR3))(n=299)
print (PDF)

pheatmap (log (hd + 1, 10),
    scale = ifelse (ncol(hd) > 2, "row", "column"),
    breaks = NA, color = palette, cluster_rows = TRUE,
    cluster_cols = ncol(hd) > 2, show_rownames = as.logical(ShowRownames),
    border_color = NA, filename = PDF, heigh = PDFheigh, width = PDFwidth)
