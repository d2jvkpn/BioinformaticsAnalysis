suppressPackageStartupMessages(library(GO.db))
library(rjson)

gc <- c("GO:0008150", "GO:0005575", "GO:0003674")
names(gc) <- c("BP", "CC", "MF")
gn <- sapply(gc, function(x) GOTERM[[x]]@Term)

anc <- c(GOBPANCESTOR, GOCCANCESTOR, GOMFANCESTOR)
names(anc) <- c("BP", "CC", "MF")

l1 <- sapply(as.list(GOTERM), function(x) return(c(x@Term, x@Ontology)))
d1 <- as.data.frame(t(as.data.frame(l1)))

GetAnc <- function(go) {
   # cat(paste0(go, "\n"))
   if (is.null(GOTERM[[go]])) { return(c(NA, NA, NA, NA)) }
   anc1 <- anc[GOTERM[[go]]@Ontology][[1]]
   jump <- c()
   a <- go
   repeat {
      s <- anc1[[a]]
      if (length(s) == 1) {break}
      a <- s[length(s)-1]
      jump <- c(a, jump)
   }

   if (length(jump) == 1) {
       jump <- c(jump, NA)
   } else if (length(jump) == 0) {
       jump <- c(go, NA)
   }

   return (c(jump[2]))
}

# d2 <- sapply(as.list(GOTERM), function(x) c(x@Term, x@Ontology, x@Definition))

for (k in names(anc)) {
    x <- sapply(as.list(anc[[k]]), function(x) paste(x, collapse=" "))
    d1[names(x), "ancestor"] <- as.character(x)
}

d1 <- cbind(rownames(d1), d1)
colnames(d1) <- c("GO_id", "term", "ontology", "ancestor")
d1$ontology <- gn[as.character(d1$ontology)]

d2 <- cbind(names(gc), as.data.frame(t(rbind(gc, gn))))
colnames(d2) <- c("abbr", "GO_id", "term")

# write.table(d2, "GO_ontology.tsv", sep="\t", quote=FALSE, row.names=FALSE)


Cell2String <- function(l) {
    if (any(is.na(l))) {return ("")}
    paste(paste(names(l), l, sep="::"), collapse=" ")
}

####
children <- c(sapply(GOBPCHILDREN, function(x) Cell2String(x)), 
    sapply(GOCCCHILDREN, function(x) Cell2String(x)),
    sapply(GOMFCHILDREN, function(x) Cell2String(x)))

d1[, "children"] <- children[rownames(d1)]

parents <- c(sapply(GOBPPARENTS, function(x) Cell2String(x)), 
    sapply(GOCCPARENTS, function(x) Cell2String(x)),
    sapply(GOMFPARENTS, function(x) Cell2String(x)))
# write.table(parents, "GO_parents.tsv", sep="\t", quote=FALSE, na="", row.names=FALSE)
d1[, "parents"] <- parents[rownames(d1)]

write.table(d1, "GO_infor.tsv", sep="\t", quote=FALSE, row.names=FALSE)
