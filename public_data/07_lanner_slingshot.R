#################Fig. 1A & Suppl.Fig 1A-B#######################################
library(SCP)
human <- readRDS("~/embryo_dataset/psd.R3.6.em.seurat.ob.rds")
human
#UpdateSeuratObject(human)
order_vec <- c(
  "Zygote", "2-4 cell", "8 cell", "Morula",
  "Prelineage", "ICM", "TE", "Early_EPI", "Early_Hypoblast",
  "EPI.PrE.INT", "Late_EPI", "Late_Hypoblast",
  "CTB", "STB", "EVT",
  "Amnion.Ecto", "Late_Amnion",
  "YSE", "AdvMes",
  "ExE_Mes", "Mesoderm", "Axial Mes", "PriS",
  "Blood_Progenitors", "Hemogenic_Endothelium",
  "Erythro-Myeloid_Progenitors", "Myeloid_Progenitors", "Erythroblasts",
  "PGC", "DE",
  "Ambiguous", "Unknown"
)
human$sub_rename_EML <- factor(human$sub_rename_EML, levels = order_vec)
DimPlot(human, reduction = "UMAP", group.by = "sub_rename_EML")
human@images <- list()
CellDimPlot(
  srt = human, group.by = c("sub_rename_EML"),
  reduction = "UMAP", theme_use = "theme_blank"
)
CellDimPlot(
  srt = human, group.by = c("devTime"),
  reduction = "UMAP", theme_use = "theme_blank"
)
#################Fig. 1B-E & Suppl.Fig 1C-D#####################################
needed <- c(
  # Hypoblast lineage
  "Early_Hypoblast", "Late_Hypoblast", "YSE", "AdvMes", "DE",
  # Hypoblast-derived mesenchyme
  "ExE_Mes",
  # Mesoderm → Hemogenic → Erythroid
  "PriS", "Mesoderm", "Axial Mes",
  "Hemogenic_Endothelium",
  "Blood_Progenitors",
  "Erythro-Myeloid_Progenitors",
  "Erythroblasts",
  #Others
  "ICM", "Early_EPI", "EPI.PrE.INT", "Late_EPI", "Amnion.Ecto", "Late_Amnion",
  "PGC", "Ambiguous"
)
Idents(human)
Idents(human) <- human@meta.data$sub_rename_EML
subset_obj <- subset(human, idents = needed)
subset_obj
DimPlot(subset_obj)
#trajectory
subset_obj <- RunSlingshot(srt = subset_obj, group.by = "sub_rename_EML", 
                           reduction = "UMAP")
CellDimPlot(subset_obj, group.by = "sub_rename_EML", reduction = "UMAP", 
            lineages = paste0("Lineage", 1:4), lineages_span = 0.1)
FeatureDimPlot(subset_obj, features = paste0("Lineage", 1:4), 
               reduction = "UMAP", theme_use = "theme_blank")
DynamicPlot(
  srt = subset_obj, lineages = c("Lineage4"), group.by = "sub_rename_EML",
  features = c("GATA6", "CDX2", "HBE1", "GYPA"),
  compare_lineages = TRUE, compare_features = FALSE
)
FeatureDimPlot(subset_obj, c("GATA6", "CDX2", "HBE1", "GYPA"), 
               reduction = "UMAP", theme_use = "theme_blank")

#################Fig. 1F-I######################################################
#Hypoblast lineage subset (Lineage 4)
subset_needed <- c(
"Early_Hypoblast", "Late_Hypoblast", "DE", "AdvMes", "ExE_Mes", "YSE"
)
hypo_obj <- subset(subset_obj, idents = needed)

hypo_obj <- NormalizeData(hypo_obj)
hypo_obj <- FindVariableFeatures(hypo_obj)

hypo_obj <- ScaleData(hypo_obj)
hypo_obj <- RunPCA(hypo_obj)
hypo_obj <- FindNeighbors(hypo_obj, dims = 1:20)
hypo_obj <- FindClusters(hypo_obj, resolution = 0.8)

hypo_obj <- RunUMAP(hypo_obj, dims = 1:20)
DimPlot(hypo_obj, group.by = "seurat_clusters", label = TRUE)

#CDX2+ Hypoblast annotation
FeaturePlot(hypo_obj, features = "CDX2", cols = c("lightgrey", "red"))
VlnPlot(hypo_obj, features = "CDX2", group.by = "seurat_clusters", pt.size = 0)
cdx2_pos <- FetchData(hypo_obj, vars = "CDX2")[,1] > 0
hypo_obj$CDX2_pos <- cdx2_pos

prop_cdx2 <- prop.table(table(hypo_obj$seurat_clusters, hypo_obj$CDX2_pos), 1)
prop_cdx2 <- as.data.frame.matrix(prop_cdx2)
prop_cdx2$cluster <- rownames(prop_cdx2)
prop_cdx2 <- prop_cdx2[order(prop_cdx2$`TRUE`, decreasing = TRUE), ]
prop_cdx2

avg_cdx2 <- AverageExpression(hypo_obj, features = "CDX2", group.by = "seurat_clusters")$RNA
avg_cdx2 <- data.frame(cluster = rownames(avg_cdx2), avg_CDX2 = avg_cdx2[,1])
avg_cdx2[order(avg_cdx2$avg_CDX2, decreasing = TRUE), ]

cdx2_cluster <- prop_cdx2$cluster[6]
cdx2_cluster

hypo_obj$Hypoblast_subtype2 <- "Other"
hypo_obj$Hypoblast_subtype2[hypo_obj$seurat_clusters == cdx2_cluster] <- "CDX2+ Hypoblast"

DimPlot(hypo_obj, group.by = "Hypoblast_subtype2", label = TRUE, pt.size = 0.5) +
  ggplot2::ggtitle(paste0("CDX2+ Hypoblast = seurat_cluster ", cdx2_cluster))

DotPlot(
  hypo_obj,
  features = c("POU5F1","SOX17","FOXA2","CXCR4","GATA6","PDGFRA","APOA1","CDX2",
               "HAND1","VIM","HBZ","HBE1","GYPA","GYPB"),
  group.by = "seurat_clusters"
) + RotatedAxis()

new.ids <- levels(hypo_obj)
new.ids[new.ids == cdx2_cluster] <- "CDX2+ Hypoblast"
hypo_obj <- RenameIdents(hypo_obj, new.ids)

CellDimPlot(
  srt = hypo_obj, group.by = Idents(hypo_obj),
  reduction = "UMAP", theme_use = "theme_blank"
)

FeatureStatPlot(hypo_obj, stat.by = c(
  "POU5F1","SOX17","FOXA2","CXCR4","GATA6","PDGFRA","APOA1","CDX2",
  "HAND1","VIM","HBZ","HBE1","GYPA","GYPB"
), group.by = Idents(hypo_obj), plot.by = "feature", stack = TRUE)

#trajectory of hypoblast subset
hypo_obj <- RunSlingshot(srt = hypo_obj, group.by = "sub_rename_EML",
                         start = "Early_Hypoblast",
                           reduction = "UMAP")
CellDimPlot(subset_obj, group.by = "sub_rename_EML", reduction = "UMAP", 
            lineages = paste0("Lineage", 1:3), lineages_span = 0.1)
FeatureDimPlot(subset_obj, features = paste0("Lineage", 1:3), 
               reduction = "UMAP", theme_use = "theme_blank")
DynamicPlot(
  srt = subset_obj, lineages = c("Lineage3"), group.by = "sub_rename_EML",
  features = c("GATA6", "CDX2", "HBE1", "GYPA"),
  compare_lineages = TRUE, compare_features = FALSE
)





