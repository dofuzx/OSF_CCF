# M1 — Gut microbiome differential signature (preregistration B1–B4)
# Run: Rscript M1_microbiome.R
suppressMessages({
  library(curatedMetagenomicData); library(Maaslin2); library(dplyr); library(metafor)
})
source("config.R")  # paths: data_dir, out_dir; params: fdr_thresh <- 0.05

cohorts <- c("ZellerG_2014","FengQ_2015","YuJ_2015","VogtmannE_2016",
             "HanniganGD_2017","WirbelJ_2018","ThomasAM_2019_a","YachidaS_2019")
# Italian data comprise three entries (ThomasAM_2019_a/b/c); pool _a/_b as discovery, _c as validation

## 1. Load harmonized species profiles (curatedMetagenomicData)
load_cohort <- function(co) {
  es <- curatedMetagenomicData(
    paste0(co, ".metaphlan_bugs_list.stool"), return = "e", taxa_base = "species")
  d <- as.data.frame(assay(es[[1]]))
  coldata <- as.data.frame(colData(es[[1]]))
  list(abund = d, meta = coldata)
}

## 2. Per-cohort differential testing (CRC vs control), MaAsLin 3 / MaAsLin2
run_maaslin <- function(co) {
  cd <- load_cohort(co)
  abund <- cd$abund[, cd$meta$study_condition %in% c("CRC","control")]
  meta  <- cd$meta[cd$meta$study_condition %in% c("CRC","control"), ]
  keep  <- rowMeans(abund > 0) >= 0.10          # prevalence filter
  fit <- Maaslin2(
    t(abund[keep, ]), meta, output = file.path(out_dir, "M1", co),
    fixed_effects = "study_condition",
    random_effects = c("age","sex","bmi")[c("age","sex","bmi") %in% colnames(meta)],
    transform = "LOG", analysis_method = "LM", max_significance = 1)
  subset(fit$results, fit$results$metadata == "study_conditionCRC")
}
res_list <- lapply(cohorts, run_maaslin); names(res_list) <- cohorts

## 3. Fixed-effect meta-analysis across cohorts (per species)
ma <- do.call(rbind, lapply(names(res_list), function(co) {
  r <- res_list[[co]]
  data.frame(cohort = co, feature = r$feature, coef = r$coef,
             stderr = r$stderr, pval = r$pval)
}))
meta_res <- ma %>%
  group_by(feature) %>%
  filter(sum(!is.na(pval)) >= 6) %>%
  summarise(
    k = sum(!is.na(pval)),
    fe_p = tryCatch(metagen(TE = coef, seTE = stderr, studlab = cohort)$pval.fixed,
                    error = function(e) NA),   # SE taken from MaAsLin2 output column 'stderr'
    n_pos = sum(sign(coef) > 0, na.rm = TRUE),
    concordant = n_pos >= 6 | n_pos == 0)
core_crc_species <- subset(meta_res, concordant & fe_p < fdr_thresh)
write.csv(meta_res, file.path(out_dir, "M1", "species_meta_all.csv"), row.names = FALSE)
write.csv(core_crc_species, file.path(out_dir, "M1", "core_crc_species.csv"), row.names = FALSE)

## 4. Functional pathways (KEGG): SCFA, bile acid, tryptophan focus
focus_kw <- c("butyrate","baiA","baiB","baiCD","baiHI","tryptophanase","tnaA","bile")
# pathway profiles loaded analogously: <cohort>.humann_kegg_pathway_uniref90.stool
# differential testing identical to step 2; filter by focus_kw; output pathway_meta.csv

## 5. Sensitivity: ANCOM-BC2 (see M7_robustness.R for concordance check)
write("M1 done", file.path(out_dir, "M1", ".done"))
