Args <- commandArgs (T)

if (length (Args) < 2 || length (Args) > 4 ) {
	print ("Arguments:  <output_prefix>  name1:input1.tsv name2:input2.tsv....")
	print ("author: d2jvkpn")
	print ("version: 0.1")
	print ("release: 2018-10-28")
	print ("project: https://github.com/d2jvkpn/BioinformaticsAnalysis")
	print ("lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html")
	q ("no", 2, FALSE)
}

prefix <- Args [1]
Sets <- strsplit(Args [2:length(Args)], ":")

dir.create ( dirname (prefix), showWarnings = FALSE, recursive = TRUE)

#### Venn data
readSet <- function (i) {
    de <- i[1]
    g <- read.delim(i[2], stringsAsFactors=F)[, 1]
    d <- data.frame(id=g, de=1)
    colnames(d)[2] <- de
    return (d)
}

Union <- Reduce(function(x, y) merge(x, y, by=1, all=TRUE), lapply(Sets, readSet))
Union[ is.na (Union) ] <- 0
Intersection <- Union[apply(Union[, 2:ncol(Union)], 1, function(x) all(x == 1)), ]

print(sprintf("Venn output: %s...", prefix))

write.table (Union, paste0(prefix, ".union.tsv"), row.names=F, quote=F, sep="\t")

write.table (Intersection, paste0(prefix, ".intersection.tsv"), row.names=F, 
quote=F, sep="\t")

#### Venn plot
if (ncol(Union) < 7) {
    library(VennDiagram)
    futile.logger::flog.threshold (futile.logger::ERROR, name = "VennDiagramLogger")

    Venn <- lapply(Union[, 2:ncol(Union)], function(x) which(x==1) )

    p <- venn.diagram (Venn, margin=0.3, lwd=0.25, 
        fill=rainbow(length(Venn)), lty=0,
		# col = rainbow(length(Venn)), print.mode=c("raw","percent"),
        cex=0.8, cat.cex=0.8, scaled=10, print.mode=c("raw"),
        category.name=names(Venn), filename= NULL)

    pdf (paste0(prefix, '.pdf'), height=8, width=8)
    grid.draw(p)
    dev.off()

} else {
    library(shape)
    library(grDevices)

    Venn <- Union[, 2:ncol(Union)]
    n <- ncol(Venn)
    cols <- rainbow(n, alpha=0.6)
    corenum <- length(which(rowSums(Venn) == n))

    pdf( paste0(prefix, '.pdf') )
    par(mar=c(1, 1, 1, 1))
    emptyplot(c(-3, 3))

    for (i in 1:n) {
        angle <- i*360/n
		Angle <- i*2*pi/n

        plotellipse(rx=1, ry=5.0/n, mid = c(cos(Angle), sin(Angle)), 
        angle=angle, col=cols[i], type='n', lwd=1)

        num <- length( which(Venn[,i]==1 & rowSums(Venn) ==1) )
        text(1.5*cos( Angle), 1.5*sin(Angle), num, cex=0.8-n/30*0.2)
        adj <- ifelse(angle > 90, ifelse(angle < 270, 1, 0), 0)

        srt <- ifelse(angle > 90, ifelse(angle < 270, angle-180, angle), angle)

        text(2*cos(Angle), 2*sin( Angle), colnames(Venn)[i], 
        adj=adj, srt=srt, cex=0.6-n/40*0.2)
    }

    plotcircle (r=0.8, col='#008FFF', type='n', lwd=1)
    text (0, 0, paste0 ('core\n', corenum), font=2)
    dev.off ()
}
