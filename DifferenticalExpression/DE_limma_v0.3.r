# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

# source('https://bioconductor.org/biocLite.R')
# options(BioC_mirror='https://mirrors.ustc.edu.cn/bioc/')
# biocLite('limma')

library(limma)
library(reshape2)

Args <- commandArgs(TRUE)

tsv <- Args[1]
prefix <- Args[2]
treated <- strsplit(Args[3], ",")[[1]]
untreated <- strsplit(Args[4], ",")[[1]]
LOG2FCMIN <- log(as.numeric(Args[5]), 2)
DESIG <- Args[6]
DESIGMIN <- as.numeric(Args[7])

d <- read.delim(tsv, stringsAsFactors=FALSE)
ed <- d[, c(untreated, treated)]

ct <- c("untreated", "treated")
group <- factor(rep(ct, each=length(treated), levels=ct))
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group","", colnames(design))

fit <- lmFit(ed, design)

contrast.matrix <- makeContrasts(treated-untreated, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

de0 <- topTable(fit2, coef=NULL, number=nrow(ed), genelist=fit$genes, 
adjust.method="fdr", sort.by="none", resort.by=NULL, p.value=1, lfc=0, 
confint=FALSE)

de <- cbind(d[,1], de0, d[, 2:ncol(d)])

colnames(de)[c(1,5,6)] <- c(colnames(d)[1], "Pvalue", "FDR")

ds <- de[abs(de$logFC) > LOG2FCMIN & de[, DESIG] < DESIGMIN, ]

write.table(de, paste0(prefix, ".DE.tsv"), sep='\t', quote=FALSE,
row.names=FALSE)

write.table(ds, paste0(prefix, '.DEsig.tsv'), sep='\t', quote=FALSE,
row.names=FALSE)
