# M2 — scRNA-seq: FOLR2+ TAM annotation, trajectory, regulons (preregistration C1–C4)
suppressMessages({library(Seurat); library(harmony); library(monocle3); library(SeuratWrappers)})
source("config.R")  # gse_dirs, seed <- 20261122

## 1. Load & QC each GSE, then integrate
## NOTE: GSE132465 and GSE178341 ship PROCESSED matrices only (raw FASTQ withheld;
## GSE178341 raw is in dbGaP phs002407). For these, load h5/CSV counts directly;
## doublet detection (scDblFinder/DoubletFinder) still runs on the count matrix.
## Not recoverable: ambient-RNA correction (CellBender needs unfiltered raw h5) and
## re-quantification with a unified Cell Ranger version — stratify these cohorts in
## sensitivity analyses (M7) and record QC differences in outputs.
qc <- function(gse_dir, processed = FALSE) {
  if (processed) {
    mat <- Read10X_h5(file.path(gse_dir, "processed.h5"))  # or read counts CSV per GSE layout
    o <- CreateSeuratObject(mat, min.cells = 3)
    o[["percent.mt"]] <- PercentageFeatureSet(o, "^MT-")
    o <- subset(o, nCount_RNA > 200 & percent.mt < 15)
    o <- scDblFinder(o)   # works on filtered count matrices
  } else {
    mat <- Read10X(gse_dir)                               # raw 10x matrix
    o <- CreateSeuratObject(mat, min.cells = 3)
    o[["percent.mt"]] <- PercentageFeatureSet(o, "^MT-")
    o <- subset(o, nCount_RNA > 200 & percent.mt < 15)
    o <- scDblFinder(o)
  }
}
so_list <- lapply(gse_dirs, qc)
so <- merge(so_list[[1]], so_list[-1]); so <- NormalizeData(so) |>
  FindVariableFeatures() |> ScaleData() |>
  RunPCA() |> RunHarmony("dataset") |>
  RunUMAP(reduction = "harmony", dims = 1:30) |>
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.6)

## 2. Myeloid subset and FOLR2+ TAM annotation (locked criteria, C1)
## Identify macrophage clusters by canonical markers, then subset
markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25)
mac_clusters <- unique(markers$cluster[markers$gene %in% c("C1QA","C1QB","CD68","APOE")])
macro <- subset(so, idents = mac_clusters)
Idents(macro) <- "seurat_clusters"
tams <- subset(macro, idents = mac_clusters)

tam_obj <- tams@assays$RNA$data
folr2_hi <- tam_obj["FOLR2", ] > quantile(tam_obj["FOLR2", ], 0.9)
res_marker <- t(tam_obj[c("MRC1","VCAM1","LYVE1","PLTP"), ]) > 1
n_pos <- rowSums(res_marker) >= 2
mono_low <- colMeans(t(tam_obj[c("S100A8","S100A9","FCN1"), ]) > 1) < 0.3
tams$FOLR2TAM <- folr2_hi & n_pos & mono_low

## 3. Signature (top 30 markers) + cross-tissue validation vs GSE192935
sig <- markers[markers$cluster %in% mac_clusters, ] |>
  subset(p_val_adj < 0.01 & avg_log2FC > 0.5) |>
  (\(m) m[order(m$avg_log2FC, decreasing = TRUE), ])()
folr2_sig <- head(unique(sig$gene[sig$gene %in% rownames(tam_obj)]), 30)
# Hypergeometric overlap vs breast FOLR2+ signature (Ramos 2022) loaded from ref file
breast_sig <- readLines(file.path(data_dir, "GSE192935", "folr2_breast_signature.txt"))
test <- phyper(length(intersect(folr2_sig, breast_sig)) - 1,
               length(breast_sig), ncol(tam_obj) - length(breast_sig),
               length(folr2_sig), lower.tail = FALSE)
write(data.frame(sig = list(folr2_sig), overlap_p = test),
      file.path(out_dir, "M2", "folr2_signature.rds"))

## 4. Trajectory (monocle3) + SCENIC regulons; RNA velocity via scVelo (external)
cds <- as.cell_data_set(tams); cds <- preprocess_cds(cds) |> reduce_dimension() |>
  cluster_cells() |> learn_graph() |> order_cells()
tams$pseudotime <- pseudotime(cds)
saveRDS(tams, file.path(out_dir, "M2", "tams_annotated.rds"))
# SCENIC (pySCENIC via reticulate) and scVelo run as separate Python steps per SOP
write("M2 done", file.path(out_dir, "M2", ".done"))
