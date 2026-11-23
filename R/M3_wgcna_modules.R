# M3 — WGCNA co-expression modules on TCGA bulk (preregistration F1)
suppressMessages({library(WGCNA); library(edgeR); library(TCGAbiolinks)})
source("config.R"); allowWGCNAThreads(); set.seed(20261122)

## 1. TCGA-COAD/READ counts -> vst-ish logCPM matrix
query <- GDCquery(project = "TCGA-COAD", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(query); GDCprepare(query) -> gdc_co   # repeat for READ; rbind
counts <- as.matrix(SummarizedExperiment::assay(gdc_co, "unfiltered_counts"))
clin   <- as.data.frame(SummarizedExperiment::colData(gdc_co))
y <- DGEList(counts); keep <- filterByExpr(y); y <- y[keep, ]
lcpm <- t(log2(edgeR::cpm(y) + 1))

## 2. Soft threshold (scale-free R^2 > 0.8)
powers <- 1:20; sft <- pickSoftThreshold(lcpm, powerVector = powers)
soft_power <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq > 0.8)]

## 3. Modules + module-trait correlation with FOLR2+ score & microbiome proxies
net <- blockwiseModules(lcpm, power = soft_power, TOMType = "unsigned",
                        minModuleSize = 30, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE)
traits <- data.frame(
  folr2_score = readRDS(file.path(out_dir, "M4", "folr2_bulk_score.rds")),
  fnuc_proxy  = readRDS(file.path(out_dir, "M4", "fnuc_proxy_score.rds")),  # F. nucleatum/dca proxy from M1
  stage = clin$ajcc_pathologic_tumor_stage)
MEs <- orderMEs(moduleEigengenes(lcpm, net$colors)$eigengenes)
mod_cor <- cor(MEs, traits, use = "pairwise"); mod_p <- corPvalueStudent(mod_cor, nrow(lcpm))

## 4. Hub genes (kME top 10) of modules significant for folr2_score
kME <- signedKME(lcpm, MEs)
write.csv(list(mod_cor = mod_cor, mod_p = mod_p, kME = kME),
          file.path(out_dir, "M3", "wgcna_results.csv"))
write("M3 done", file.path(out_dir, "M3", ".done"))
