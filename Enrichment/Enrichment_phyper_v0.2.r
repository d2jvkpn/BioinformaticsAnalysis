# Project: https://github.com/d2jvkpn/BioinformaticsAnalysis

Args <- commandArgs(TRUE)

target <- read.delim(Args[1], stringsAsFactors=FALSE)
map2 <- read.delim(Args[2], stringsAsFactors=FALSE)
n <- as.integer(Args[3])
out <- Args[4]

if (n<2) {
    print(sprintf("invalid column number %d", n))
    q(save = "default", 1)
}

tids <- unique(target[target[,1]!="",1])
map2 <- map2[map2[,1] != "" & map2[,n] != "", ]
map2 <- map2[!duplicated(map2),]

td <- map2[, n:ncol(map2)]
colnames(td)[1] <- "term"
td <- td[!duplicated(td$term), ]

if( length(tids[tids %in% map2[,1]]) == 0) {
   print(sprintf("%s no match term in %s found", Args[1], Args[2]))
   q(save = "default", 0)
}

dtt <- map2[map2[, 1] %in% tids, ]

d <- aggregate(dtt[,1], by=list(dtt[, n]), function(x) length(x))

colnames(d) <- c("term", "q")
d$m <- sapply(d[,1], function(x) nrow(map2[map2[,n] == x,]))
d$n <- length(unique(map2[,1])) - d$m
d$k <- length(unique(dtt[,1]))

d$count <- apply(d[, 2:5], 1, function(x) paste(x[1], x[2], x[3], x[4],
    sep=", "))

d$score <- (d$q/d$k)/(d$m/(d$m + d$n))
d$Pvalue <- apply(d[, 2:5], 1, function(x) phyper(x[1], x[2], x[3], x[4]))
d$FDR <- p.adjust(d$Pvalue, method = "fdr")

d$ids <- aggregate(dtt[,1], by=list(dtt[, n]),
    function(x) paste(x, collapse="; "))[[2]]

d1 <- merge(td[td$term %in% d$term, ],
    d[, c("term", "count","score","Pvalue","FDR", "ids")], 
    by.x="term", by.y="term")

d1 <- d1 [order(d1$Pvalue), ]
colnames(d1)[1] <- colnames(map2)[n]

dir.create(dirname(d1), showWarnings = FALSE, recursive = TRUE)

write.table(d1, out, sep="\t", quote=FALSE, row.names=FALSE)
print(sprintf("saved %s", out))
