

# Python

# import pandas as pd
# df = pd.read_csv("../data/GSE136106_Dumdie_countdata.txt", sep="\t", index_col=0)
# n_samples = 28
# raw = df.iloc[:, :28]
# norm = df.iloc[:, 28:56]

# norm.columns = norm.columns.str.replace(".1", "", regex=False)

# icm_cols = [c for c in raw.columns if "_M" in c]

# raw_icm = raw[icm_cols]
# norm_icm = norm[icm_cols]

# meta = pd.DataFrame({
#     "sample": raw_icm.columns
# })
# meta["group"] = meta["sample"].apply(
#     lambda x: "HP" if x.startswith("G_") else "LP"
# )

# # gene symbol to gene name
# import mygene
# mg = mygene.MyGeneInfo()
# gene_ids = raw_icm.index.tolist()
# res = mg.querymany(
#     gene_ids,
#     scopes="ensembl.gene",
#     fields="symbol",
#     species="human"
# )
# mapping = pd.DataFrame(res)
# mapping = mapping[["query", "symbol"]].dropna()
# mapping_dict = dict(zip(mapping["query"], mapping["symbol"]))

# raw_icm["symbol"] = raw_icm.index.map(mapping_dict).copy()
# norm_icm["symbol"] = norm_icm.index.map(mapping_dict).copy()

# raw_icm = raw_icm.dropna(subset=["symbol"])
# norm_icm = norm_icm.dropna(subset=["symbol"])

# raw_icm = raw_icm.set_index("symbol")
# norm_icm = norm_icm.set_index("symbol")



# R
setwd('/mnt/yiming/nfs_share/blood/notebooks/')

library(DESeq2)

counts <- read.csv("../data/raw_icm_unique.csv", row.names=1)
coldata <- read.csv("../data/coldata.csv", row.names=1)

all(colnames(counts) == rownames(coldata))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$padj),]

write.csv(as.data.frame(resOrdered),
          "../data/DESeq2_results.csv")

summary(res)


library(ggplot2)
library(ggrepel) 
library(dplyr)

plot_volcano_deseq2 <- function(de, 
                                title = "HP vs LP Volcano Plot", 
                                lfc_threshold = 0.5, 
                                pval_threshold = 0.05, 
                                use_fdr = TRUE, 
                                highlight_genes = NULL) {
  
  df <- as.data.frame(de)
  
  if (use_fdr) {
    df$p_use <- df$padj
    p_label <- "FDR"
  } else {
    df$p_use <- df$pvalue
    p_label <- "p-value"
  }
  
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$p_use), ]
  df$neg_log10 <- -log10(df$p_use + 1e-300)
  
  df$group <- "Not significant"
  df$group[df$p_use <= pval_threshold & df$log2FoldChange >= lfc_threshold] <- "Up in HP"
  df$group[df$p_use <= pval_threshold & df$log2FoldChange <= -lfc_threshold] <- "Up in LP"
  
  df$group <- factor(df$group, levels = c("Up in HP", "Up in LP", "Not significant"))
  
  up_hp_genes <- rownames(df[df$group == "Up in HP", ])
  up_lp_genes <- rownames(df[df$group == "Up in LP", ])
  
  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10, color = group)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up in HP" = "red", 
                                  "Up in LP" = "blue", 
                                  "Not significant" = "lightgray")) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
    labs(title = title,
         x = "log2 Fold Change (HP - LP)",
         y = paste0("-log10(", p_label, ")"),
         color = "Significance") +
    theme_minimal() +
    theme(legend.position = "right")
  
  if (!is.null(highlight_genes)) {
    genes_to_plot <- df[rownames(df) %in% highlight_genes, ]
    
    if (nrow(genes_to_plot) > 0) {
      p <- p + 
        geom_point(data = genes_to_plot, color = "black", size = 3, shape = 8) +
        geom_text_repel(data = genes_to_plot, aes(label = rownames(genes_to_plot)),
                        color = "black", size = 3, fontface = "bold")
    }
  }
  
  cat("Total genes:", nrow(df), "\n")
  cat("Significant:", sum(df$group != "Not significant"), "\n")
  cat("Up in HP:", length(up_hp_genes), "\n")
  cat("Up in LP:", length(up_lp_genes), "\n")
  
  print(p)
  
  return(list(up_hp_genes = up_hp_genes, up_lp_genes = up_lp_genes))
}

results_list <- plot_volcano_deseq2(resOrdered, lfc_threshold = 1.0)

my_genes <- c('CITED4')
plot_volcano_deseq2(res, highlight_genes = my_genes)

