# OSF_CCF — 结直肠癌 TDO 与巨噬细胞多组学研究资源

> 本仓库托管 OSF_CCF 项目的**脱敏研究数据、分析流程与文档**，并通过 Zenodo 归档获取 DOI。

## 项目概述

本项目围绕结直肠癌（CRC）肿瘤来源类器官（TDO）与肿瘤相关巨噬细胞（TAM），整合公开多组学队列（肠道宏基因组、单细胞转录组、bulk RNA-seq），开展微生物组—免疫微环境—FOLR2⁺ 巨噬细胞特征的系统分析。分析流程涵盖数据下载、单细胞注释、WGCNA 模块、机器学习特征构建、机制中介分析、分子对接准备与稳健性验证。

**所有数据均来自公开数据库或经脱敏处理，不含可直接识别患者身份的信息。**

## 目录结构

```
OSF_CCF/
├── Document/
│   ├── Study Design Proposal.pdf            研究设计与总体方案
│   ├── Data Download SOP.pdf                数据下载标准操作流程
│   ├── FOLR2 Pipeline Specification.pdf     FOLR2 分析流程规范
│   ├── Wet-Lab Validation Framework.pdf     湿实验验证框架
│   ├── data_inventory.csv                   数据源清单（队列、登录号、样本量）
│   └── spatial_manual_checklist_template.csv 空间转录组人工核查模板
├── R/
│   ├── M1_microbiome.R                       宏基因组数据处理
│   ├── M2_scrna_annotation.R                 单细胞注释
│   ├── M3_wgcna_modules.R                    WGCNA 模块分析
│   ├── M4_ml_signature.R                     机器学习特征构建
│   ├── M5_mechanism_mediation.R              机制与中介分析
│   ├── M6_docking_prep.R                     分子对接前处理
│   ├── M7_robustness.R                       稳健性验证
│   └── spatial_download_backup.R             空间数据下载备份脚本
├── deviation/                               偏离记录（每次偏离一个 Markdown 文件）
├── README.md
├── CITATION.cff
└── .gitignore
```

## 数据来源

数据源清单见 `Document/data_inventory.csv`，涵盖 ENA / SRA / DDBJ（宏基因组）、GEO（单细胞与 bulk）、TCGA / GDC 等公开库。每条记录含队列、国家、登录号、样本量与分析角色。使用第三方数据请遵守各原始数据库的授权与署名要求。

## 版本与 DOI

- https://doi.org/10.5281/zenodo.24459024

## 引用

若使用本资源，请引用 `CITATION.cff` 中对应的版本 DOI。

## 数据使用与伦理声明

- 本仓库内容仅用于**学术研究**，不得用于临床诊断或商业用途（除非获得明确授权）。
- 数据来自公开数据库或经脱敏处理，仓库不包含可直接识别个人身份的信息（PHI）。
- 使用者应遵守原始数据提供方的伦理审批、知情同意与数据使用条款。
- 若发现任何潜在敏感信息泄露，请联系维护者。

## 许可

- 分析脚本（`/R`）：建议 MIT License（请按机构要求确认后添加 `LICENSE` 文件）。
- 文档与数据：建议 CC BY 4.0（请按机构与数据来源要求确认）。

## 维护者

`[Fu Zexian、Hebei University of Engineering、fuzexian@hebeu.edu.cn、ORCID：https://orcid.org/0000-0002-6945-0486]`
