Args <- commandArgs(T)

if (length(Args) == 0 || Args[1] == "-h" || Args[1] == "--help") {
	print("Arguments:  <expression.tsv>  <outputPrefix>  <threshold>  <title>  <group.tsv>")
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
group <- Args[5]

d <- read.delim (tsv, check.names=F, quote="", row.names=1)
gp <- read.delim(group, row.names=1)

if (FALSE %in% c(sort(colnames(d)) == sort(rownames(gp))) ) {
	print("Error: sample names doesn't match the group.")
	q ()
}

library(ggplot2)
library(reshape2)
library (pheatmap)
library (gridExtra)

m <- 1:length(unique (gp[,1]))
names(m) <- unique (gp[,1])
d <- d[, order (m[gp[colnames(d), 1]])]

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

# dsum ggplot2

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
md1$sample <- as.character(md1$sample)
md1$group <- gp[md1$sample, 1]


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

bplot <- ggplot(md1, aes(x=sample, y=log10title, fill=group)) + 
	coord_cartesian (y = c(log(threshold, 10), max(md$log10title))) + 
	geom_boxplot (notch=FALSE, notchwidth = 0.3, alpha=0.8) +
	ggtitle (paste (title, "boxplot")) + 
	xlab ("") +
	ylab(bquote ("log" [10] ~ .(title))) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5),
#	legend.position="none",
	axis.text.x=element_text(angle = 45, hjust=1))

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
	aes(x=sample, y=log10title, fill=group))+ 
	coord_cartesian (y = c(log(threshold, 10), max(md$log10title))) + 
	geom_violin() + 
	xlab("") +
	# geom_jitter(size=0.5, alpha=0.6) +
	ggtitle(paste(title, "violinplot")) + 
	ylab(bquote("log" [10] ~ .(title))) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5),
#	legend.position="none",
	axis.text.x=element_text(angle = 45, hjust=1)) +
	stat_summary(fun.data=data_summary)

ggsave (paste0(prefix, ".violinplot.pdf"), plot=vplot)


####
print("Densityplot...")
pdf (paste0(prefix, ".densityplot.pdf"), onefile=TRUE)

for (g in unique(gp[,1])) {
	s <- colnames(d) [g == gp[,1]]
	if (length(s) < 2) { next }

	dplot <- ggplot(md1[md1$sample %in% s, ], 
	aes(x=log10title, colour=sample)) +
	#coord_cartesian(x = c(-3, max(md$log10name))) +
	geom_density() + 
	ggtitle (paste (g, title, "density")) +
	xlab(bquote("log" [10] ~ .(title))) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5) )

	print (dplot)
}

dev.off()


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

####
print("Correlation scatter plot...")

CorScatter <- function(d, p, threshold, title) {
	p <- as.character(p)
	d1 <- d[apply(d[, p], 1, min) >= threshold, p]
	d1[d1==0] <- 1E-6
	d1 <- log(d1, 10)
	colnames(d1) <- c('a', 'b')

	pearson <- round (cor(d1[,1], d1[,2],  method = "pearson"), 3)
	kendall <- round (cor(d1[,1], d1[,2],  method = "kendall"), 3)
	spearman <- round (cor(d1[,1], d1[,2],  method = "spearman"), 3)

	ann <- paste (paste0('pearson = ', pearson),
	paste0('kendall = ', kendall),
	paste0('spearman = ', spearman), sep='\n')

	pl <- ggplot (d1, aes (a, b)) + 
	geom_point (colour="blue", alpha=0.2, size=0.5) +
	ggtitle(paste(p[2],"~", p[1], "scatterplot" )) +
	xlab (bquote('log' [10] ~ .(title) ~ ", " ~ .(p[1]))) + 
	ylab (bquote('log' [10] ~ .(title) ~ ", " ~ .(p[2]))) +
	# annotate ("text", -Inf, Inf, label = ann, hjust = 0, vjust = 1) +
	annotate ("text", min(d1[,1]), max(d1[,2]), 
	label = ann, hjust = 0, vjust = 1) +
	theme_bw() +
	theme (aspect.ratio=1, plot.title = element_text(hjust = 0.5))

	return (pl)
}

for (g in unique(gp[,1])) {
	s <- colnames(d) [g == gp[,1]]
	if(length(s) ==1) { next }

	pdf (paste0 (dirname(prefix), "/CorScatter_", g, ".pdf"), onefile = TRUE)

	for (p in as.data.frame(combn(s,2))) {
		print (CorScatter(d, as.character(p), threshold, title))
	}

	dev.off()
}
