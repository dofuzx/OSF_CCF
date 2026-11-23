# OSF_CCF — Multi-omics Research Resources for Colorectal Cancer TDO and Macrophages

> This repository hosts the **de-identified research data, analysis workflows, and documentation** of the OSF_CCF project, archived on Zenodo with a DOI.

## Overview

This project focuses on colorectal cancer (CRC) tumor-derived organoids (TDO) and tumor-associated macrophages (TAM). It integrates public multi-omics cohorts (gut metagenomes, single-cell transcriptomes, and bulk RNA-seq) to systematically dissect the microbiome – immune microenvironment – FOLR2⁺ macrophage signature axis. The analysis workflows cover data download, single-cell annotation, WGCNA module analysis, machine-learning signature construction, mechanistic mediation analysis, molecular docking preparation, and robustness validation.

**All data are either obtained from public databases or de-identified, and contain no directly identifiable patient information.**

## Directory Structure

```
OSF_CCF/
├── Document/
│   ├── Study Design Proposal.pdf              Study design and overall protocol
│   ├── Data Download SOP.pdf                  Standard operating procedure for data download
│   ├── FOLR2 Pipeline Specification.pdf       FOLR2 analysis pipeline specification
│   ├── Wet-Lab Validation Framework.pdf       Wet-lab validation framework
│   ├── data_inventory.csv                     Data source inventory (cohorts, accessions, sample sizes)
│   └── spatial_manual_checklist_template.csv  Spatial transcriptomics manual QC checklist template
├── R/
│   ├── M1_microbiome.R                        Metagenomic data processing
│   ├── M2_scrna_annotation.R                  Single-cell annotation
│   ├── M3_wgcna_modules.R                     WGCNA module analysis
│   ├── M4_ml_signature.R                      Machine-learning signature construction
│   ├── M5_mechanism_mediation.R               Mechanistic and mediation analysis
│   ├── M6_docking_prep.R                      Molecular docking preparation
│   ├── M7_robustness.R                        Robustness validation
│   └── spatial_download_backup.R              Spatial data download backup script
├── deviation/                                 Deviation records (one Markdown file per deviation)
├── README.md                                  Chinese version
├── README_en.md                               English version
├── CITATION.cff
└── .gitignore
```

## Data Sources

The data inventory is provided in `Document/data_inventory.csv`, covering public repositories such as ENA / SRA / DDBJ (metagenomes), GEO (single-cell and bulk RNA-seq), and TCGA / GDC. Each entry records the cohort, country, accession number, sample size, and analytical role. When using third-party data, please comply with the licensing and attribution requirements of the original databases.

## Versions and DOI

- https://doi.org/10.5281/zenodo.24459742

## Citation

If you use this resource, please cite the corresponding version DOI in `CITATION.cff`.

## Data Use and Ethics Statement

- The content of this repository is for **academic research only** and must not be used for clinical diagnosis or commercial purposes (unless explicitly authorized).
- The data are obtained from public databases or de-identified; this repository contains no protected health information (PHI) that could directly identify an individual.
- Users must comply with the ethical approvals, informed consent, and data-use terms of the original data providers.
- If any potential sensitive information leakage is identified, please contact the maintainer.

## License

- Analysis scripts (`/R`): MIT License recommended (please add a `LICENSE` file after confirming with your institution).
- Documentation and data: CC BY 4.0 recommended (please confirm against institutional requirements and the terms of the original data sources).

## Maintainers

`[Fu Zexian, Hebei University of Engineering, fuzexian@hebeu.edu.cn, ORCID: https://orcid.org/0000-0002-6945-0486]`
