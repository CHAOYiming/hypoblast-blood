#devtools::install_github("liulab-dfci/MAGeCKFlute")
library(MAGeCKFlute)

gdata = ReadRRA("difftest.gene_summary.txt")
head(gdata)

sdata = ReadsgRRA("difftest.sgrna_summary.txt")
head(sdata)

countsummary = read.delim("crispr.countsummary.txt", check.names = FALSE)

BarView(countsummary,
        x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads")

BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")

MapRatesView(countsummary)

#Fig. 6C & F
gdata$Rank = rank(gdata$Score)
ScatterView(gdata, x = "Rank", y = "Score", label = "id",
            top = 5, auto_cut_y = TRUE, ylab = "Log2FC",
            groups = c("top", "bottom"))
