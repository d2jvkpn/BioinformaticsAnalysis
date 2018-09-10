Args <- commandArgs(T)

if (length(Args) == 0 || Args[1] == "-h" || Args[1] == "--help") {
	print("Arguments:  <expression.tsv>  <outputPrefix>  <threshold>  <title>")
	print("author: d2jvkpn")
	print("version: 0.3")
	print("release: 2018-09-10")
	print("project: https://github.com/d2jvkpn/BioinformaticsAnalysis")
	print("lisense: GPLv3 (https://www.gnu.org/licenses/gpl-3.0.en.html")
	q()
}

tsv <- Args[1]
prefix <- Args[2]
threshold <- as.double (Args[3])
title <- Args[4]

d <- read.delim (tsv, check.names=F, quote="", row.names=1)

library(ggplot2)
library(reshape2)
library (pheatmap)
library (gridExtra)

dis <- function (x, threshold) {
	x <- x[x >= threshold]
	s = quantile (x, probs = c(0, 0.25, 0.5, 0.75, 1.0))
	s["count"] <- length(x)
	s["mean"] <- mean(x)
	s["std"] <- sd(x)
	s["threshold"] <- threshold

	s <- s[c ("threshold", "count", "mean", "std", "0%", "25%", "50%", "75%", 
	"100%")]

	names (s)[5:9] <- c("min", "Q1", "median", "Q3", "max")

	return (round (s, 3))
}

dir.create (dirname (prefix), showWarnings=FALSE, recursive=TRUE)

print("Summary...")

top10 <- t (apply(d, 2, function(x) sort(x, decreasing=TRUE)[1:10]))

n <- apply(d, 2, function(x) order(x, decreasing=TRUE)[1:10])
top10gene <- t (apply (n, 2, function(x) rownames(d)[x]))

colnames (top10gene) <- c("top1", "top2", "top3", "top4", "top5", "top6", 
"top7", "top8", "top9", "top10")

for (i in 1:nrow(top10gene) ) {
	for (j in 1:ncol(top10gene) ) {
		top10gene[i,j] <- paste0(top10gene[i,j] , " (", top10[i, j], ")")
	}
}

write.table (top10gene, paste0(prefix, ".top10_genes.tsv"), col.names=NA, 
sep="\t", quote=F)

dsum <- t (apply(d, 2, function (x) dis(x, 0.01)))
write.table (dsum, paste0(prefix, ".summary.tsv"), col.names=NA, sep="\t", quote=F)


####
md <- melt (d)
colnames (md) <- c("sample", "title")
md$log10title <- log(md$title, 10)

md[which(md[,2] < 0.01), "interval"] <- "[0, 0.01)"
md[which(md[,2] >= 0.01 & md[,2] < 0.1), "interval"] <- "[0.01, 0.1)"
md[which(md[,2] >= 0.1 & md[,2] < 1), "interval"] <- "[0.1, 1)"
md[which(md[,2] >= 1 & md[,2] < 10), "interval"] <- "[1, 10)"
md[which(md[,2] >= 10 & md[,2] < 20), "interval"] <- "[10, 20)"
md[which(md[,2] >= 20 & md[,2] < 100), "interval"] <- "[20, 100)"
md[which(md[,2] >= 100), "interval"] <- "[100, Inf)"

L1 <- c("[100, Inf)", "[20, 100)", "[10, 20)", "[1, 10)", "[0.1, 1)",
"[0.01, 0.1)", "[0, 0.01)")

md[, "interval"] <- factor (md[, "interval"], levels=L1)

md1 <- md[which(md$log10title >= log(threshold, 10)), ]


###
print("Interval plot...")

interval <- table (md[, c("sample", "interval")])

write.table (interval[, ncol(interval):1],
paste0(prefix, ".interval.tsv"), col.names=NA, sep="\t", quote=F)

iplot <- ggplot(md, aes(sample, fill=interval)) +
	geom_histogram (stat="count") + 
	ggtitle(paste(title, "interval")) + 
	xlab("") +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5),
	axis.text.x=element_text(angle = 45, hjust=1))

ggsave (paste0(prefix, ".interval.pdf"), plot=iplot)


####
print("Boxplot...")

bplot <- ggplot(md1, aes(x=sample, y=log10title, fill=sample)) + 
	coord_cartesian (y = c(log(threshold, 10), max(md$log10title))) + 
	geom_boxplot (notch=FALSE, notchwidth = 0.3, alpha=0.8) +
	ggtitle (paste (title, "boxplot")) + 
	xlab ("") +
	ylab(bquote ("log" [10] ~ .(title))) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5),
	axis.text.x=element_text(angle = 45, hjust=1), 
		legend.position="none") +
	guides(fill=guide_legend(title=""))

ggsave(paste0(prefix, ".boxplot.pdf"), plot=bplot)


####
print("Violinplot...")

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

vplot <- ggplot(md1,
	aes(x=sample, y=log10title, fill=sample))+ 
	coord_cartesian (y = c(log(threshold, 10), max(md$log10title))) + 
	geom_violin() + 
	xlab("") +
	# geom_jitter(size=0.5, alpha=0.6) +
	ggtitle(paste(title, "violinplot")) + 
	ylab(bquote("log" [10] ~ .(title))) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5),
	axis.text.x=element_text(angle = 45, hjust=1), 
		legend.position="none") +
	guides(fill=guide_legend(title="")) +
	stat_summary(fun.data=data_summary)

ggsave (paste0(prefix, ".violinplot.pdf"), plot=vplot)


####
print("Densityplot...")

dplot <- ggplot(md1, aes(x=log10title, colour=sample)) +
	#coord_cartesian(x = c(-3, max(md$log10name))) +
	geom_density() + 
	ggtitle (paste(title, "density")) +
	xlab(bquote("log" [10] ~ .(title))) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5) ) +
	guides(fill=guide_legend(title=""))

ggsave(paste0(prefix, ".densityplot.pdf"), plot=dplot)


####
if ((ncol(d)) < 2) {q()}

print("Correlation heatmap...")

CorrM <- function (d, threshold, method) {
	k <- ncol(d)
	corrm <- data.frame (matrix(ncol=k, nrow=k), row.names = colnames(d))

	colnames(corrm) <- colnames(d)

	for(i in 1:(ncol(corrm)-1)) {

		for(j in (i+1):ncol(corrm)) {
			x <- d[,c(i,j)]
			x <- x[apply(x, 1, min) >= threshold, ]
			x[x==0] <- 1E-6
			y <- cor(log(x[,1], 10), log(x[,2], 10), method = method)
			corrm[i, j] <- y
			corrm[j, i] <- y
		}
	}

	corrm[is.na(corrm)] <- 1

	return (corrm)
}

for (m in c("pearson", "kendall", "spearman")) {
	corrm <- CorrM(d, threshold, m)

	pheatmap (corrm, cluster_rows=FALSE, 
	cluster_cols=FALSE, display_numbers=TRUE,
	main=paste("Correlation in", m), number_format="%.3f",
	file=paste0 (prefix, ".corheatmap_", m, ".pdf"))
}
