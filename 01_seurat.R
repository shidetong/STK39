# ============================================================
# Part 1. 
# ============================================================

setwd('/work/home/20260401_STK39')

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

if (!requireNamespace("harmony", quietly = TRUE)) {
  install.packages("harmony")
}
library(harmony)

# -------------------- 1.  --------------------
data_dir <- '/work/home/stk39/data'
sample_folders <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)

if (length(sample_folders) < 25) {
  warning(paste("仅", length(sample_folders), "个样本文件夹"))
}

if (length(sample_folders) > 25) {
  sample_folders <- sample_folders[1:25]
}

sample_names <- basename(sample_folders)
seurat_list <- list()

for (i in seq_along(sample_folders)) {
  tryCatch({
    sample_data <- Read10X(sample_folders[i])

    if (is.null(sample_data) || ncol(sample_data) == 0) {
      next
    }

    seurat_obj <- CreateSeuratObject(
      counts = sample_data,
      project = sample_names[i],
      min.cells = 3,
      min.features = 200
    )

    if (ncol(seurat_obj) == 0) {
      next
    }

    seurat_obj$sample <- sample_names[i]
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
      seurat_obj,
      pattern = "^MT-"
    )

    seurat_list[[sample_names[i]]] <- seurat_obj
  }, error = function(e) {
    warning(
      paste(
        "读取样本", sample_names[i], "出错:",
        e$message
      )
    )
  })
}

if (length(seurat_list) == 0) {
  stop("没有成功读取任何样本数据")
}

# -------------------- 2. --------------------
age_meta_file <- '/work/home/stk39/data/sample_metadata.csv'

if (file.exists(age_meta_file)) {
  sample_meta <- read.csv(
    age_meta_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  colnames(sample_meta) <- c('sample_name', 'age')

  for (sample_name in names(seurat_list)) {
    if (sample_name %in% sample_meta$sample_name) {
      age <- sample_meta$age[
        sample_meta$sample_name == sample_name
      ]
      seurat_list[[sample_name]]$age <- age
    } else {
      seurat_list[[sample_name]]$age <- "Unknown"
    }
  }
} else {
  for (sample_name in names(seurat_list)) {
    seurat_list[[sample_name]]$age <- "Unknown"
  }
}

# -------------------- 3. --------------------
# 200 < nFeature_RNA < 5500
# 800 < nCount_RNA < 17000
# percent.mt < 25

qc_summary <- data.frame()

for (sample_name in names(seurat_list)) {
  seurat_obj <- seurat_list[[sample_name]]
  cells_before <- ncol(seurat_obj)

  seurat_obj <- subset(
    seurat_obj,
    subset =
      nFeature_RNA > 200 &
      nFeature_RNA < 5500 &
      nCount_RNA > 800 &
      nCount_RNA < 17000 &
      percent.mt < 25
  )

  cells_after <- ncol(seurat_obj)
  seurat_list[[sample_name]] <- seurat_obj

  qc_summary <- rbind(
    qc_summary,
    data.frame(
      sample = sample_name,
      cells_before = cells_before,
      cells_after = cells_after,
      cells_removed = cells_before - cells_after
    )
  )
}

write.csv(
  qc_summary,
  "TESTIS_QC_SUMMARY.csv",
  row.names = FALSE
)

# -------------------- 4. --------------------
integrated <- merge(
  seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)]
)

rm(seurat_list)
gc()

# -------------------- 5.  --------------------
integrated <- NormalizeData(
  integrated,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

integrated <- FindVariableFeatures(
  integrated,
  selection.method = "vst",
  nfeatures = 5000
)

integrated <- ScaleData(
  integrated,
  features = rownames(integrated)
)

# -------------------- 6.  --------------------
integrated <- RunPCA(
  integrated,
  features = VariableFeatures(integrated),
  npcs = 50,
  verbose = FALSE
)

pdf("TESTIS_ELBOW_PLOT.pdf", width = 8, height = 6)
print(
  ElbowPlot(integrated, ndims = 50) +
    ggtitle("Elbow Plot for PCA Dimension Selection")
)
dev.off()

pcs_use <- 1:30

# -------------------- 7.  --------------------
integrated <- RunHarmony(
  object = integrated,
  group.by.vars = "sample",
  reduction = "pca",
  dims.use = pcs_use,
  assay.use = "RNA",
  reduction.save = "harmony",
  verbose = FALSE
)

# -------------------- 8.  --------------------
integrated <- RunUMAP(
  integrated,
  reduction = "harmony",
  dims = pcs_use,
  verbose = FALSE
)

integrated <- FindNeighbors(
  integrated,
  reduction = "harmony",
  dims = pcs_use
)

integrated <- FindClusters(
  integrated,
  resolution = 0.4
)

saveRDS(
  integrated,
  "TESTIS_SEURAT_CLUSTERED.rds"
)

# -------------------- 9.  --------------------
# Cluster UMAP
p_cluster <- DimPlot(
  integrated,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
) +
  ggtitle("Testis Cells Clustering (Harmony Batch Correction)") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  )

ggsave(
  "TESTIS_UMAP_CLUSTERS.pdf",
  p_cluster,
  width = 10,
  height = 8
)

# Sample UMAP：
p_batch <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "sample",
  pt.size = 0.5
) +
  ggtitle("Testis Cells by Sample") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  )

ggsave(
  "TESTIS_UMAP_BATCH.pdf",
  p_batch,
  width = 10,
  height = 8
)

# -------------------- 10.  --------------------
cluster_summary <- data.frame(
  Cluster = names(table(integrated$seurat_clusters)),
  Cell_Count = as.numeric(
    table(integrated$seurat_clusters)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_summary,
  "TESTIS_CLUSTER_SUMMARY.csv",
  row.names = FALSE
)
