# M5 — Mechanism: cross-cohort association, mediation, cell-cell comm, survival (D, E1, E3)
suppressMessages({library(Maaslin2); library(survival); library(CellChat); library(nichenetr); library(mediation)})
source("config.R"); set.seed(20261122)

## 1. D1: pathway/species abundance vs FOLR2+ TAM score (cohort-level, mixed model)
dat <- readRDS(file.path(out_dir, "M4", "association_table.rds"))  # rows=cohorts/stages
fit_mixed <- lmerTest::lmer(folr2_score ~ pathway_dca + (1 | cohort), data = dat)
spearman <- cor.test(dat$pathway_dca, dat$folr2_score, method = "spearman")

## 2. D2: paired atlas (Roelands 2023) tissue microbiome x FOLR2+ score
paired <- readRDS(file.path(data_dir, "Roelands2023", "paired_tme_microbiome.rds"))
paired_cor <- apply(paired$taxa, 2, function(x)
  cor.test(x, paired$folr2_score, method = "spearman")$p.value)

## 3. D3: TME correlates + CIBERSORTx (via immunedeconv / Docker)
folr2_score_bulk <- readRDS(file.path(out_dir, "M4", "folr2_bulk_score.rds"))
expr <- readRDS(file.path(data_dir, "TCGA", "tpm.rds"))
# CIBERSORTx absolute mode (Docker) -> cell fractions; correlate CD8/Treg vs folr2_score

## 4. E1: mediation panel -> FOLR2+ score -> OS (pre-specified panel, ranked by proportion mediated)
med_data <- readRDS(file.path(out_dir, "M4", "mediation_input.rds"))
mediators <- c("pathway_dca", "fnuc_abundance", "pathway_scfas", "pathway_tryp", "fnev_proxy")
med_out_list <- lapply(mediators, function(x) {
  mod_m <- lm(folr2_score ~ .(x) + age + sex + stage, data = med_data)   # pseudo: construct per-mediator formula
  mod_y <- coxph(Surv(OS.time, OS) ~ folr2_score + .(x) + age + sex + stage, data = med_data)
  mediation::mediate(mod_m, mod_y, treat = x, mediator = "folr2_score", boot = TRUE, sims = 1000)
})
# additionally: restricted cubic splines for SCFA non-linearity (rms::rcs);
# negative-control taxa (non-CRC oral commensals) tested for specificity

## 5. E3: KM (median split) + Cox with MSI/TMB/IPS interactions
dd <- med_data; dd$high <- as.numeric(dd$folr2_score > median(dd$folr2_score))
km_fit <- survfit(Surv(OS.time, OS) ~ high, data = dd); survdiff(Surv(OS.time, OS) ~ high, data = dd)
cox_full <- coxph(Surv(OS.time, OS) ~ high + age + sex + stage + MSI +
                    high:MSI, data = dd)

## 6. CellChat/NicheNet on scRNA TAM <-> T cell pairs (per dataset, aggregated)
# CellChat: subset epithelial/myeloid/T; overexpressed LR pairs FOLR2+ TAM -> CD8/Treg
# NicheNet: prior-based ligand activity for LGALS9, TGF-beta, IL-10 programs
saveRDS(list(mixed = fit_mixed, med = med_out_list, cox = cox_full),
        file.path(out_dir, "M5", "mechanism_results.rds"))
write("M5 done", file.path(out_dir, "M5", ".done"))
