library(SCP)
library(MouseGastrulationData)

data("mouse_gastrulation", package = "MouseGastrulationData")
mouse <- mouse_gastrulation

meta_cols <- colnames(mouse@meta.data)
stage_col <- meta_cols[meta_cols %in% c("main_celltype", "stage", "Stage", "timepoint", 
                                        "Timepoint", "embryo_stage", "embryonic_day", "day")]
if (length(stage_col) == 0) {
  message("No obvious stage column found. Skipping stage subset; using all cells.")
} else {
  stage_col <- stage_col[1]
  keep <- grepl("^E6\\.?5|^E7\\.?0|^E7\\.?5", as.character(mouse@meta.data[[stage_col]]))
  if (any(keep)) obj <- subset(mouse, cells = rownames(mouse@meta.data)[keep])
}

DefaultAssay(mouse) <- "RNA"
mouse <- NormalizeData(mouse)
mouse <- FindVariableFeatures(mouse, nfeatures = 3000)
mouse <- ScaleData(mouse, features = VariableFeatures(obj))
mouse <- RunPCA(mouse, npcs = 30, verbose = FALSE)
mouse <- FindNeighbors(mouse, dims = 1:30)
mouse <- FindClusters(mouse, resolution = 0.6)
mouse <- RunUMAP(mouse, dims = 1:30)

#################Suppl.Fig 1E-G#################################################
CellDimPlot(
  srt = mouse, group.by = c("main_celltype"),
  reduction = "UMAP", theme_use = "theme_blank"
)

markers <- c("Sox17","Gata6","Gata4","Cdx2","Vim",
             "Apoa2","Gypa","Gata1","Hba-a1","Hbb-bs")
pF <- DotPlot(
  mouse,
  features = markers,
  group.by = ct_col
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  ggtitle("Marker gene expression level") +
  theme(plot.title = element_text(hjust = 0.5))

mouse <- RunSlingshot(srt = mouse, group.by = "main_celltype",
                         reduction = "umap")
FeatureDimPlot(mouse, features = paste0("Lineage", 1:5), 
               reduction = "umap", theme_use = "theme_blank")


