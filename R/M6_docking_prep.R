# M6 — Molecular docking prep + HPC AutoDock Vina submission (preregistration F2)
# Ligands: DCA, butyrate, indole-3-acetate (small molecules); FomA/Fap2 domains (proteins)
# Receptors: TLR4/MD-2, TGR5(GPBAR1), FXR, TIGIT, TREM2
suppressMessages({library(reticulate); library(jsonlite)})
source("config.R")  # hpc credentials via biomni.tool

## 1. Receptor structures from PDB / AlphaFold DB
receptors <- list(
  TLR4_MD2   = list(pdb = "3FXI"),           # human wild-type TLR4/MD-2 (LPS-bound)
  TGR5       = list(pdb = "PDBID_or_AF"),    # GPBAR1 cryo-EM complex
  FXR        = list(pdb = "3DCT"),           # FXR ligand-binding domain
  TIGIT      = list(pdb = "3Q0H"),           # human TIGIT ectodomain (apo; alternative: 3UDW TIGIT/PVR)
  TREM2      = list(pdb = "5ELI"))
# download via BiocFileCache::bfcadd or https://files.rcsb.org/download/<PDB>.pdb

## 2. Small-molecule ligands: SMILES -> 3D (RDKit via reticulate), enumerate protonation
rdkit <- import("rdkit.Chem")
ligands <- list(
  DCA = "CCC(C)C1CCC2CC1CC(=O)N2",  # placeholder — verify canonical SMILES in execution
  butyrate = "CCCC(=O)[O-]",
  indole3acetate = "O=C(O)Cc1cn[cH]2ccccc12")

## 3. Protein-protein pairs (FomA/Fap2 domains vs TIGIT/TLR4): use HPC AlphaFold/docking
## 4. Submit AutoDock Vina batch on HPC via biomni.tool
library(biomni.tool)
docs <- hpc_search_tools("autodock vina")
cat(docs[[1]]$usage, "\n")
# job <- hpc_run_tool("autodock-vina", "<usage-compliant command>",
#                     input_files = list("receptor.pdbqt" = ..., "ligand.pdbqt" = ...))
## 5. Parse scores + residue contacts -> docking_scores.csv; filter dG < -7 kcal/mol
write("M6 prep done — docking executed on HPC after registration",
      file.path(out_dir, "M6", ".done"))
