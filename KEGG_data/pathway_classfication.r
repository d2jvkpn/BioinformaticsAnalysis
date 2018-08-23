Args <- commandArgs(TRUE)

d = read.delim(Args[1], stringsAsFactors=FALSE)
PDF = Args[2]

colnames(d)[1:2] <- c("L2", "category")

d[d$category=="Environmental Information Processing", 
"category"] <- "Environmental Information\nProcessing"

d[d$category=="Genetic Information Processing", 
"category"] <- "Genetic Information\nProcessing"

d$L2 <- factor(d$L2, levels=rev(d$L2))
d$category <- factor(d$category, levels=unique(d$category))

library(ggplot2)

p <- ggplot(d, aes(x=L2, y=gene_count, fill = category)) +
    geom_bar(stat="identity") +
    ggtitle("Pathway Classification") +
    theme_bw(base_size=10) +
    theme( plot.title = element_text(hjust = 0.5, size=10), 
        legend.position="bottom",
        plot.margin= unit(c(0.5, 2, 0.5, 0.5), "cm")) +
    scale_fill_discrete("") +
    scale_y_continuous(limits=c(0, max(d$gene_count) * 1.2) ) +
    geom_text(aes(label=paste0 (d$gene_count, " (", d$pathway_count, ")")), 
        hjust=-0.1, size=2.5) +
    xlab("") +
    coord_flip() +
    labs(y = "gene_count (pathway_count)")

ggsave(PDF, p, height=8, width=8)
