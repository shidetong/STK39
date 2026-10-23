# ============================================================
# Part 2. 
# ============================================================

setwd('/work/home/20260401_STK39')

library(Seurat)
library(dplyr)
library(ggplot2)

# -------------------- 1.  --------------------
integrated <- readRDS(
  "TESTIS_SEURAT_CLUSTERED.rds"
)

DefaultAssay(integrated) <- "RNA"

integrated <- JoinLayers(integrated)

Idents(integrated) <- integrated$seurat_clusters

# -------------------- 2.  --------------------
all_markers <- FindAllMarkers(
  integrated,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

if (nrow(all_markers) == 0) {
  all_markers <- FindAllMarkers(
    integrated,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.1
  )
}

if (nrow(all_markers) == 0) {
  stop("未找到cluster marker")
}

write.table(
  all_markers,
  "TESTIS_CLUSTER_MARKERS.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  )

write.table(
  top10_markers,
  "TESTIS_TOP10_MARKERS.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------- 3.  --------------------
testis_marker_genes <- c(
  "PIWIL4", "TSPAN33", "ID4", "TCF3",
  "UTF1", "MKI67", "KIT", "HMGN2", "STRA8",
  "SYCP3", "SPO11", "MEIOB", "TOP2A",
  "SYCP1", "SPAG6", "SPATA8", "CETN3",
  "ZPBP", "ACRV1", "CCDC168", "SPACA1",
  "FAM186A", "PRM3", "TNP2", "TSSK6",
  "SPATA3", "LELP1", "FATE1", "SOX9",
  "WT1", "AMH", "ACTA2", "MYH11",
  "WFDC1", "PTCH1", "IGF1", "DLK1",
  "NEGR1", "NOTCH3", "RGS5", "VWF",
  "CD34", "CD14", "C1QC", "CD3D",
  "CD69", "CD79A", "IGHM", "JCHAIN"
)

available_markers <- testis_marker_genes[
  testis_marker_genes %in% rownames(integrated)
]

write.table(
  available_markers,
  "TESTIS_AVAILABLE_MARKERS.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# -------------------- 4.  --------------------
pdf(
  "TESTIS_MARKER_DOTPLOT.pdf",
  width = 14,
  height = 6
)

p_marker <- DotPlot(
  integrated,
  features = unique(available_markers),
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 10)
  ) +
  ggtitle("Testis Marker Gene Dot Plot")

print(p_marker)
dev.off()

# -------------------- 5. --------------------
cluster_to_celltype <- list(
  '0'  = "Endothelial cell",
  '1'  = "TPC",
  '2'  = "Leydig cell",
  '3'  = "Round spermatids",
  '4'  = "Round spermatids",
  '5'  = "SMC",
  '6'  = "Macrophage",
  '7'  = "Endothelial cell",
  '8'  = "Early primary spermatocytes",
  '9'  = "Elongating spermatids",
  '10' = "Elongating spermatids",
  '11' = "Elongating spermatids",
  '12' = "Leydig cell",
  '13' = "Undifferentiated spermatogonia",
  '14' = "Late primary spermatocytes",
  '15' = "Differentiating spermatogonia",
  '16' = "Lymphocyte",
  '17' = "Sertoli cell",
  '18' = "Endothelial cell",
  '19' = "Macrophage",
  '20' = "Endothelial cell"
)

all_clusters <- as.character(
  unique(integrated$seurat_clusters)
)

unannotated <- setdiff(
  all_clusters,
  names(cluster_to_celltype)
)

if (length(unannotated) > 0) {
  for (clust in unannotated) {
    cluster_to_celltype[[clust]] <- paste0(
      "Unknown_",
      clust
    )
  }
}

current_clusters <- as.character(
  integrated$seurat_clusters
)

new_celltypes <- sapply(
  current_clusters,
  function(x) cluster_to_celltype[[x]]
)

names(new_celltypes) <- colnames(integrated)

integrated <- AddMetaData(
  integrated,
  metadata = new_celltypes,
  col.name = "cell_type"
)

Idents(integrated) <- integrated$cell_type

# -------------------- 6.  --------------------
if ("age" %in% colnames(integrated@meta.data)) {
  integrated$age <- as.numeric(
    as.character(integrated$age)
  )

  integrated$age_group <- sapply(
    integrated$age,
    function(x) {
      if (is.na(x)) return("Unknown")
      if (x >= 20 & x < 30) "20s"
      else if (x >= 30 & x < 40) "30s"
      else if (x >= 40 & x < 50) "40s"
      else if (x >= 50 & x < 60) "50s"
      else if (x >= 60 & x < 70) "60s"
      else "other"
    }
  )
} else {
  integrated$age_group <- "Unknown"
  integrated$age <- NA
}

# -------------------- 7.  --------------------
p_annotation <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.3
) +
  ggtitle("Testis Cell Type Annotation") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  )

ggsave(
  "TESTIS_MANUAL_ANNOTATION_UMAP.pdf",
  p_annotation,
  width = 12,
  height = 9
)

# -------------------- 8.  --------------------
cell_type_counts <- table(
  integrated$cell_type
)

cell_type_distribution <- data.frame(
  Cell_Type = names(cell_type_counts),
  Cell_Count = as.numeric(cell_type_counts),
  stringsAsFactors = FALSE
)

write.csv(
  cell_type_distribution,
  "TESTIS_CELLTYPE_DISTRIBUTION.csv",
  row.names = FALSE
)

saveRDS(
  integrated,
  "TESTIS_FINAL_ANNOTATED.rds"
)
