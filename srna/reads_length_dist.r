library(ggplot2)

# fs <- c(Sys.glob("blastn/*.length.tsv"), Sys.glob("clean/*.length.tsv"))
Args <- commandArgs(TRUE)
fs = Args

for (f in fs) {
  d <- read.delim(f)
  PDF <- sub(".tsv", ".pdf", f)

  p <- ggplot(d, aes(x=length, y=count)) + 
    geom_bar(stat="identity", fill="#4682B4") +
    ## theme_bw() +
    theme(plot.title = element_text(hjust = 0.5) ) +
    ggtitle("Length distribution") + 
    ylab("") + xlab("") + scale_x_continuous(breaks = seq(min(d$length), max(d$length), by = 2)) +
    geom_text(aes(label=d$count), vjust=-0.5, size=2)

  ggsave(PDF, p)
}
