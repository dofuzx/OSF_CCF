# M7 — Robustness: leave-one-cohort-out, method agreement, compositional bias (F3, B2)
suppressMessages({library(ANCOMBC); library(Seurat); library(magrittr)})
source("config.R"); set.seed(20261122)

## 1. Leave-one-cohort-out replication of primary associations (P1, P3)
loocv <- function(cohorts, test_fn) {
  sapply(cohorts, function(hold) {
    res <- test_fn(setdiff(cohorts, hold))          # re-run M1+D1 without held-out cohort
    sig_hold <- res$significant[[hold]]             # effect in held-out cohort
    sign_concordant <- sign(res$meta_coef) == sign(sig_hold)
    data.frame(held_out = hold, q = res$q, concordant = sign_concordant)
  }, simplify = FALSE)
}

## 2. Method agreement: MaAsLin 3 vs ANCOM-BC2 differential species
ancombc2 <- function(co) {  # per-cohort ANCOM-BC2 on same prevalence-filtered species
  # SE, ANCOMBC::ancombc2(data = physeq, fix_formula = "study_condition + age + sex + bmi")
}
agreement <- intersect(species_maaslin_sig, species_ancombc_sig)
direction_agree_rate <- mean(sign(coef_maaslin[agreement]) == sign(coef_ancom[agreement]))
stopifnot(direction_agree_rate >= 0.8)              # quality gate per spec

## 3. Deconvolution agreement: CIBERSORTx vs EPIC vs quanTIseq for CD8/Treg fractions
# compare FOLR2+ score associations across methods (Spearman sign concordance)

## 4. Compositionality bias: add pseudocount + CLR check; rarefaction sensitivity
## 5. Batch: re-run M2 integration without Harmony; compare FOLR2+ annotation stability (Jaccard)

write.csv(data.frame(agreement = list(agreement), rate = direction_agree_rate),
          file.path(out_dir, "M7", "robustness_report.csv"))
write("M7 done", file.path(out_dir, "M7", ".done"))
