# merge BM and eAML data for eAML1

library(dplyr)
library(ggplot2)
library(Seurat)

AML44.BM <- readRDS("./data/scRNA_seq/objects/AML44_BM.rds")
AML44.eAML <- readRDS("./data/scRNA_seq/objects/20230117_eAML1.rds")

AML44 <- merge(AML44.BM, AML44.eAML)
AML44 <- JoinLayers(AML44)

source("./scRNA_seq/R/00_read_clonotypes.R")

read_clonotypes("./data/scRNA_seq/20230207_TCR_eAML1_with_marrow.xlsx", AML44, TCR.clones.file = "./data/scRNA_seq/objects/20250429_eAML1.combined.TCR.clones.csv")

AML44 <- NormalizeData(AML44)
AML44 <- ScaleData(AML44)

saveRDS(AML44, file = "./data/scRNA_seq/objects/20250429_eAML1.combined.rds")
