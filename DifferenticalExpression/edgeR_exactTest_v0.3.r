# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

library("edgeR")

Args <- commandArgs(T)

# Args <- c("GSE92322/GSE92322.infor1.tsv", "T:T1,T2,T3,T4,T5",
# "N:N1,N2,N3,N4,N5", "GSE92322/GSE92322")
tsv <- Args[1]
treated <- Args[2]
untreated <- Args[3]
prefix <- Args[4]

counts <- read.delim(tsv, check.names=F)

s1 <- strsplit(strsplit(treated, ":")[[1]][2], ",")[[1]]
g1 <- strsplit(treated, ":")[[1]][1]
s2 <- strsplit(strsplit(untreated, ":")[[1]][2], ",")[[1]]
g2 <- strsplit(untreated, ":")[[1]][1]
samples <- c(list(s2), list(s1))
names(samples) <- c(g2, g1)

group <- factor(rep(names(samples), times=sapply(samples, length)),
    levels=names(samples))

if (basename(prefix) == "" ) {
    prefix = sprintf("%s/%s_vs_%s", prefix, g1, g2)
} else {
    prefix = sprintf("%s.%s_vs_%s", prefix, g1, g2)
}

y <- DGEList(counts = counts[, unlist(samples)], group = group, genes=counts[,1])
y <- calcNormFactors(y)

keep <- rowSums( cpm(y) > 1 ) >= 2
y <- y[keep, keep.lib.sizes=FALSE]

if (length(unlist(samples)) > 2) {
   y <- estimateDisp(y)
   fit <- glmFit(y)
   et <- topTags(exactTest(y), n=Inf, adjust.method="fdr", sort.by="none")
} else {
   bcv <- 0.05

   et <- topTags(exactTest(y, dispersion=bcv*2), n=Inf, adjust.method="fdr",
       sort.by="none")
}

ed <- et$table

colnames(ed)[1:4] <- c(colnames(counts)[1], "log2FoldChange", "log2CPM", 
    "Pvalue")

ed$regulated <- sapply(ed$log2FoldChange, function(x) ifelse(x>0, "up", "down"))

dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)

tsv <- paste0(prefix, ".edgeR.tsv")

write.table(ed, tsv, sep="\t", row.names=FALSE, na="", quote=F)
    print(sprintf("Saved %s", tsv))
