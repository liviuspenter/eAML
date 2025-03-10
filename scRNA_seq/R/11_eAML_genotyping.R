# analyze nanoranger data

library(ggplot2)
library(Seurat)
# need to get from github https://github.com/liviuspenter/nanoranger.R
library(nanoranger.R)
source("./scRNA_seq/R/00_eAML.colors.R")

eAML.data <- readRDS("./data/objects/20241021_eAML_combined.scvi.rds")

# read individual pileups

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.1_pileup_NPM1.csv.gz"))
eAML1.1.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.1_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML1.1.NPM1$bc <- paste0("eAML1.1_", eAML1.1.NPM1$bc, "-1")
eAML1.1.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML1.1.NPM1$bc]
write.csv2(eAML1.1.NPM1, file = "./data/scRNA_seq/genotyping/eAML1.1_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.1_pileup_STAG2.csv.gz"))
eAML1.1.STAG2 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.1_pileup_STAG2.csv.gz", REF = "G", CONSENSUS = -1, FILTER = 10)
eAML1.1.STAG2$bc <- paste0("eAML1.1_", eAML1.1.STAG2$bc, "-1")
eAML1.1.STAG2$manual.cluster <- eAML.data$manual.cluster[eAML1.1.STAG2$bc]
write.csv2(eAML1.1.STAG2, file = "./data/scRNA_seq/genotyping/eAML1.1_STAG2.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.2_pileup_NPM1.csv.gz"))
eAML1.2.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.2_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML1.2.NPM1$bc <- paste0("eAML1.2_", eAML1.2.NPM1$bc, "-1")
eAML1.2.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML1.2.NPM1$bc]
write.csv2(eAML1.2.NPM1, file = "./data/scRNA_seq/genotyping/eAML1.2_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.2_pileup_STAG2.csv.gz"))
eAML1.2.STAG2 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.2_pileup_STAG2.csv.gz", REF = "G", CONSENSUS = -1, FILTER = 10)
eAML1.2.STAG2$bc <- paste0("eAML1.2_", eAML1.2.STAG2$bc, "-1")
eAML1.2.STAG2$manual.cluster <- eAML.data$manual.cluster[eAML1.2.STAG2$bc]
write.csv2(eAML1.2.STAG2, file = "./data/scRNA_seq/genotyping/eAML1.2_STAG2.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.3_pileup_NPM1.csv.gz"))
eAML1.3.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.3_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML1.3.NPM1$bc <- paste0("eAML1.3_", eAML1.3.NPM1$bc, "-1")
eAML1.3.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML1.3.NPM1$bc]
write.csv2(eAML1.3.NPM1, file = "./data/scRNA_seq/genotyping/eAML1.3_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.3_pileup_STAG2.csv.gz"))
eAML1.3.STAG2 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.3_pileup_STAG2.csv.gz", REF = "G", CONSENSUS = -1, FILTER = 10)
eAML1.3.STAG2$bc <- paste0("eAML1.3_", eAML1.3.STAG2$bc, "-1")
eAML1.3.STAG2$manual.cluster <- eAML.data$manual.cluster[eAML1.3.STAG2$bc]
write.csv2(eAML1.3.STAG2, file = "./data/scRNA_seq/genotyping/eAML1.3_STAG2.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.4_pileup_NPM1.csv.gz"))
eAML1.4.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.4_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML1.4.NPM1$bc <- paste0("eAML1.4_", eAML1.4.NPM1$bc, "-1")
eAML1.4.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML1.4.NPM1$bc]
write.csv2(eAML1.4.NPM1, file = "./data/scRNA_seq/genotyping/eAML1.4_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.4_pileup_STAG2.csv.gz"))
eAML1.4.STAG2 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.4_pileup_STAG2.csv.gz", REF = "G", CONSENSUS = -1, FILTER = 10)
eAML1.4.STAG2$bc <- paste0("eAML1.4_", eAML1.4.STAG2$bc, "-1")
eAML1.4.STAG2$manual.cluster <- eAML.data$manual.cluster[eAML1.4.STAG2$bc]
write.csv2(eAML1.4.STAG2, file = "./data/scRNA_seq/genotyping/eAML1.4_STAG2.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.5_pileup_NPM1.csv.gz"))
eAML1.5.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.5_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML1.5.NPM1$bc <- paste0("eAML1.5_", eAML1.5.NPM1$bc, "-1")
eAML1.5.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML1.5.NPM1$bc]
write.csv2(eAML1.5.NPM1, file = "./data/scRNA_seq/genotyping/eAML1.5_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.5_pileup_STAG2.csv.gz"))
eAML1.5.STAG2 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.5_pileup_STAG2.csv.gz", REF = "G", CONSENSUS = -1, FILTER = 10)
eAML1.5.STAG2$bc <- paste0("eAML1.5_", eAML1.5.STAG2$bc, "-1")
eAML1.5.STAG2$manual.cluster <- eAML.data$manual.cluster[eAML1.5.STAG2$bc]
write.csv2(eAML1.5.STAG2, file = "./data/scRNA_seq/genotyping/eAML1.5_STAG2.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML1.6_pileup_NPM1.csv.gz"))
eAML1.6.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML1.6_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML1.6.NPM1$bc <- paste0("eAML1.6_", eAML1.6.NPM1$bc, "-1")
eAML1.6.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML1.6.NPM1$bc]
write.csv2(eAML1.6.NPM1, file = "./data/scRNA_seq/genotyping/eAML1.6_NPM1.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML2_pileup_NPM1.csv.gz"))
eAML2.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML2_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 1)
eAML2.NPM1$bc <- paste0("eAML3.1_", eAML2.NPM1$bc, "-1")
eAML2.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML2.NPM1$bc]
write.csv2(eAML2.NPM1, file = "./data/scRNA_seq/genotyping/eAML2_NPM1.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML2_pileup_DNMT3A_783.csv.gz"))
eAML2.DNMT3A <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML2_pileup_DNMT3A_783.csv.gz", REF = "T", ALT = "C", FILTER = 50)
eAML2.DNMT3A$bc <- paste0("eAML2_", eAML2.DNMT3A$bc, "-1")
eAML2.DNMT3A$manual.cluster <- eAML.data$manual.cluster[eAML2.DNMT3A$bc]
write.csv2(eAML2.DNMT3A, file = "./data/scRNA_seq/genotyping/eAML2_DNMT3A_783.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML2_pileup_DNMT3A_844.csv.gz"))
eAML2.DNMT3A <- nanoranger.R::extract_indel("./data/scRNA_seq/genotyping/eAML2_pileup_DNMT3A_844.csv.gz", REF = "T", FILTER = 10)
eAML2.DNMT3A$bc <- paste0("eAML2_", eAML2.DNMT3A$bc, "-1")
eAML2.DNMT3A$manual.cluster <- eAML.data$manual.cluster[eAML2.DNMT3A$bc]
write.csv2(eAML2.DNMT3A, file = "./data/scRNA_seq/genotyping/eAML2_DNMT3A_783.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.1_pileup_NPM1.csv.gz"))
eAML3.1.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML3.1_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML3.1.NPM1$bc <- paste0("eAML3.1_", eAML3.1.NPM1$bc, "-1")
eAML3.1.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML3.1.NPM1$bc]
write.csv2(eAML3.1.NPM1, file = "./data/scRNA_seq/genotyping/eAML3.1_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.1_pileup_DNMT3A.csv.gz"))
eAML3.1.DNMT3A <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML3.1_pileup_DNMT3A.csv.gz", REF = "G", ALT = "A", FILTER = 10)
eAML3.1.DNMT3A$bc <- paste0("eAML3.1_", eAML3.1.DNMT3A$bc, "-1")
eAML3.1.DNMT3A$manual.cluster <- eAML.data$manual.cluster[eAML3.1.DNMT3A$bc]
write.csv2(eAML3.1.DNMT3A, file = "./data/scRNA_seq/genotyping/eAML3.1_DNMT3A.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.1_pileup_IDH1.csv.gz"))
eAML3.1.IDH1 <- nanoranger.R::extract_mutation(BC.data.file = "./data/scRNA_seq/genotyping/eAML3.1_pileup_IDH1.csv.gz", REF = "C", ALT = "T", FILTER = 20)
eAML3.1.IDH1$bc <- paste0("eAML3.1_", eAML3.1.IDH1$bc, "-1")
eAML3.1.IDH1$manual.cluster <- eAML.data$manual.cluster[eAML3.1.IDH1$bc]
write.csv2(eAML3.1.IDH1, file = "./data/scRNA_seq/genotyping/eAML3.1_IDH1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.1_pileup_KRAS.csv.gz"))
eAML3.1.KRAS <- nanoranger.R::extract_mutation(BC.data.file = "./data/scRNA_seq/genotyping/eAML3.1_pileup_KRAS.csv.gz", REF = "A", ALT = "C", FILTER = 20)
eAML3.1.KRAS$bc <- paste0("eAML3.1_", eAML3.1.KRAS$bc, "-1")
eAML3.1.KRAS$manual.cluster <- eAML.data$manual.cluster[eAML3.1.KRAS$bc]
write.csv2(eAML3.1.KRAS, file = "./data/scRNA_seq/genotyping/eAML3.1_KRAS.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.2_pileup_NPM1.csv.gz"))
eAML3.2.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML3.2_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML3.2.NPM1$bc <- paste0("eAML3.2_", eAML3.2.NPM1$bc, "-1")
eAML3.2.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML3.2.NPM1$bc]
write.csv2(eAML3.2.NPM1, file = "./data/scRNA_seq/genotyping/eAML3.2_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.2_pileup_DNMT3A.csv.gz"))
eAML3.2.DNMT3A <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML3.2_pileup_DNMT3A.csv.gz", REF = "G", ALT = "A", FILTER = 10)
eAML3.2.DNMT3A$bc <- paste0("eAML3.2_", eAML3.2.DNMT3A$bc, "-1")
eAML3.2.DNMT3A$manual.cluster <- eAML.data$manual.cluster[eAML3.2.DNMT3A$bc]
write.csv2(eAML3.2.DNMT3A, file = "./data/scRNA_seq/genotyping/eAML3.2_DNMT3A.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.2_pileup_IDH1.csv.gz"))
eAML3.2.IDH1 <- nanoranger.R::extract_mutation(BC.data.file = "./data/scRNA_seq/genotyping/eAML3.2_pileup_IDH1.csv.gz", REF = "C", ALT = "T", FILTER = 20)
eAML3.2.IDH1$bc <- paste0("eAML3.2_", eAML3.2.IDH1$bc, "-1")
eAML3.2.IDH1$manual.cluster <- eAML.data$manual.cluster[eAML3.2.IDH1$bc]
write.csv2(eAML3.2.IDH1, file = "./data/scRNA_seq/genotyping/eAML3.2_IDH1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML3.2_pileup_KRAS.csv.gz"))
eAML3.2.KRAS <- nanoranger.R::extract_mutation(BC.data.file = "./data/scRNA_seq/genotyping/eAML3.2_pileup_KRAS.csv.gz", REF = "A", ALT = "C", FILTER = 10)
eAML3.2.KRAS$bc <- paste0("eAML3.2_", eAML3.2.KRAS$bc, "-1")
eAML3.2.KRAS$manual.cluster <- eAML.data$manual.cluster[eAML3.2.KRAS$bc]
write.csv2(eAML3.2.KRAS, file = "./data/scRNA_seq/genotyping/eAML3.2_KRAS.csv", quote = F)



nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.1_pileup_NPM1.csv.gz"))
eAML5.1.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML5.1_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML5.1.NPM1$bc <- paste0("eAML5.1_", eAML5.1.NPM1$bc, "-1")
eAML5.1.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML5.1.NPM1$bc]
write.csv2(eAML5.1.NPM1, file = "./data/scRNA_seq/genotyping/eAML5.1_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.1_pileup_DNMT3A.csv.gz"))
eAML5.1.DNMT3A <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML5.1_pileup_DNMT3A.csv.gz", REF = "G", ALT = "A", FILTER = 10)
eAML5.1.DNMT3A$bc <- paste0("eAML5.1_", eAML5.1.DNMT3A$bc, "-1")
eAML5.1.DNMT3A$manual.cluster <- eAML.data$manual.cluster[eAML5.1.DNMT3A$bc]
write.csv2(eAML5.1.DNMT3A, file = "./data/scRNA_seq/genotyping/eAML5.1_DNMT3A.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.1_pileup_JAK2.csv.gz"))
eAML5.1.JAK2 <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML5.1_pileup_JAK2.csv.gz", REF = "G", ALT = "T", FILTER = 10)
eAML5.1.JAK2$bc <- paste0("eAML5.1_", eAML5.1.JAK2$bc, "-1")
eAML5.1.JAK2$manual.cluster <- eAML.data$manual.cluster[eAML5.1.JAK2$bc]
write.csv2(eAML5.1.JAK2, file = "./data/scRNA_seq/genotyping/eAML5.1_JAK2.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.1_pileup_IDH1.csv.gz"))
eAML5.1.IDH1 <- nanoranger.R::extract_mutation(BC.data.file = "./data/scRNA_seq/genotyping/eAML5.1_pileup_IDH1.csv.gz", REF = "C", ALT = "T", FILTER = 20)
eAML5.1.IDH1$bc <- paste0("eAML5.1_", eAML5.1.IDH1$bc, "-1")
eAML5.1.IDH1$manual.cluster <- eAML.data$manual.cluster[eAML5.1.IDH1$bc]
write.csv2(eAML5.1.IDH1, file = "./data/scRNA_seq/genotyping/eAML5.1_IDH1.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.2_pileup_NPM1.csv.gz"))
eAML5.2.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML5.2_pileup_NPM1.csv.gz", REF = "C", CONSENSUS = 4, FILTER = 20)
eAML5.2.NPM1$bc <- paste0("eAML5.2_", eAML5.2.NPM1$bc, "-1")
eAML5.2.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML5.2.NPM1$bc]
write.csv2(eAML5.2.NPM1, file = "./data/scRNA_seq/genotyping/eAML5.2_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.2_pileup_DNMT3A.csv.gz"))
eAML5.2.DNMT3A <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML5.2_pileup_DNMT3A.csv.gz", REF = "G", ALT = "A", FILTER = 10)
eAML5.2.DNMT3A$bc <- paste0("eAML5.2_", eAML5.2.DNMT3A$bc, "-1")
eAML5.2.DNMT3A$manual.cluster <- eAML.data$manual.cluster[eAML5.2.DNMT3A$bc]
write.csv2(eAML5.2.DNMT3A, file = "./data/scRNA_seq/genotyping/eAML5.2_DNMT3A.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.2_pileup_JAK2.csv.gz"))
eAML5.2.JAK2 <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML5.2_pileup_JAK2.csv.gz", REF = "G", ALT = "T", FILTER = 10)
eAML5.2.JAK2$bc <- paste0("eAML5.2_", eAML5.2.JAK2$bc, "-1")
eAML5.2.JAK2$manual.cluster <- eAML.data$manual.cluster[eAML5.2.JAK2$bc]
write.csv2(eAML5.2.JAK2, file = "./data/scRNA_seq/genotyping/eAML5.2_JAK2.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML5.2_pileup_IDH1.csv.gz"))
eAML5.2.IDH1 <- nanoranger.R::extract_mutation(BC.data.file = "./data/scRNA_seq/genotyping/eAML5.2_pileup_IDH1.csv.gz", REF = "C", ALT = "T", FILTER = 20)
eAML5.2.IDH1$bc <- paste0("eAML5.2_", eAML5.2.IDH1$bc, "-1")
eAML5.2.IDH1$manual.cluster <- eAML.data$manual.cluster[eAML5.2.IDH1$bc]
write.csv2(eAML5.2.IDH1, file = "./data/scRNA_seq/genotyping/eAML5.2_IDH1.csv", quote = F)


nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML6_pileup_NPM1.csv.gz"))
eAML6.NPM1 <- nanoranger.R::extract_indel(BC.data.file = "./data/scRNA_seq/genotyping/eAML6_pileup_NPM1.csv.gz", REF = "T", CONSENSUS = 4, FILTER = 20)
eAML6.NPM1$bc <- paste0("eAML6_", eAML6.NPM1$bc, "-1")
eAML6.NPM1$manual.cluster <- eAML.data$manual.cluster[eAML6.NPM1$bc]
write.csv2(eAML6.NPM1, file = "./data/scRNA_seq/genotyping/eAML6_NPM1.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML6_pileup_NPM1_VUS.csv.gz"))
eAML6.NPM1_VUS <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML6_pileup_NPM1_VUS.csv.gz", REF = "A", ALT = "G", FILTER = 10)
eAML6.NPM1_VUS$bc <- paste0("eAML6_", eAML6.NPM1_VUS$bc, "-1")
eAML6.NPM1_VUS$manual.cluster <- eAML.data$manual.cluster[eAML6.NPM1_VUS$bc]
write.csv2(eAML6.NPM1_VUS, file = "./data/scRNA_seq/genotyping/eAML6_NPM1_VUS.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML6_pileup_DNMT3A.csv.gz"))
eAML6.DNMT3A <- nanoranger.R::extract_mutation("./data/scRNA_seq/genotyping/eAML6_pileup_DNMT3A.csv.gz", REF = "C", ALT = "T", FILTER = 10)
eAML6.DNMT3A$bc <- paste0("eAML6_", eAML6.DNMT3A$bc, "-1")
eAML6.DNMT3A$manual.cluster <- eAML.data$manual.cluster[eAML6.DNMT3A$bc]
write.csv2(eAML6.DNMT3A, file = "./data/scRNA_seq/genotyping/eAML6_DNMT3A.csv", quote = F)

nanoranger.R::knee_plot(data.table::fread("./data/scRNA_seq/genotyping/eAML6_pileup_IDH1.csv.gz"))
eAML6.IDH1 <- nanoranger.R::extract_mutation(BC.data.file = "./data/scRNA_seq/genotyping/eAML6_pileup_IDH1.csv.gz", REF = "G", ALT = "A", FILTER = 20)
eAML6.IDH1$bc <- paste0("eAML6_", eAML6.IDH1$bc, "-1")
eAML6.IDH1$manual.cluster <- eAML.data$manual.cluster[eAML6.IDH1$bc]
write.csv2(eAML6.IDH1, file = "./data/scRNA_seq/genotyping/eAML6_IDH1.csv", quote = F)

combined.mutations <- data.frame()
for (f in list.files("./data/scRNA_seq/genotyping/", pattern = "*.csv$", full.names = T)) {
  boo <- as.data.frame(read.table(file = f, sep = ";", header = T))
  boo$sample <- stringr::str_split_fixed(gsub(basename(f), pattern = ".csv", replacement = ""), pattern = "_", n = 2)[, 1]
  boo$target <- stringr::str_split_fixed(gsub(basename(f), pattern = ".csv", replacement = ""), pattern = "_", n = 2)[, 2]
  combined.mutations <- rbind(combined.mutations, boo)
}

combined.mutations$manual.cluster <- eAML.data$manual.cluster[combined.mutations$bc]

targets <- unique(combined.mutations$target)
combined.mutations <- combined.mutations %>%
  dplyr::select(-c(X, alt, ref, vaf)) %>%
  group_by(bc) %>%
  tidyr::pivot_wider(names_from = "target", values_from = "mutated")
combined.mutations$mutated <- apply(combined.mutations, 1, FUN = function(x) {
  if (length(which(x[targets] == "mutated")) > length(which(x[targets] == "wildtype"))) {
    "mutated"
  } else {
    "wildtype"
  }
})
combined.mutations <- as.data.frame(combined.mutations)
rownames(combined.mutations) <- combined.mutations$bc
write.csv(file = "./data/scRNA_seq/20241022_eAML_genotyping.csv", combined.mutations)

# combined.mutations = read.csv2('./data/20241022_eAML_genotyping.csv', sep = ',') %>% as.data.frame()
combined.mutations$sample <- stringr::str_split_fixed(rownames(combined.mutations), pattern = "_", n = 2)[, 1]

# plot mutational burden across celltypes
boo <- combined.mutations %>%
  group_by(sample, manual.cluster, mutated) %>%
  tally() %>%
  filter(!is.na(manual.cluster)) %>%
  tidyr::pivot_wider(names_from = "mutated", values_from = "n") %>%
  mutate(mutated = tidyr::replace_na(mutated, 0)) %>%
  mutate(wildtype = tidyr::replace_na(wildtype, 0)) %>%
  mutate(mut.freq = mutated / (mutated + wildtype))

boo <- boo[-which((boo$wildtype + boo$mutated < 3)), ]

boo$manual.cluster <- factor(boo$manual.cluster, levels = c(
  "endothelium", "fibroblast", "keratinocyte", "melanocyte", "smooth muscle",
  "myeloid", "macrophage", "T", "NK", "B cell", "erythroid", "megakaryocytic"
))
ggplot(data = boo[which(boo$manual.cluster != "erythroid"), ], aes(x = manual.cluster, y = 100 * mut.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, linewidth = 0.5, size = 0.5, color = "black") +
  geom_point(aes(color = manual.cluster), size = 0.5) +
  scale_y_continuous("% mutated") +
  scale_color_manual(values = eAML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20241022_eAML_combined_mutated.svg", width = 2.5, height = 2)

ggplot(data = boo[which(boo$manual.cluster != "erythroid"), ], aes(x = manual.cluster, y = 100 * mut.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, linewidth = 0.5, size = 0.5, color = "black") +
  geom_point(aes(color = manual.cluster), size = 1) +
  scale_y_continuous("% mutated") +
  scale_color_manual(values = eAML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250227_eAML_combined_mutated.svg", width = 3, height = 2)


eAML.data$mutated <- combined.mutations[colnames(eAML.data), "mutated"]
p <- DimPlot(eAML.data, group.by = "mutated", cols = c("wildtype" = "blue", "mutated" = "firebrick"), na.value = "grey90", raster = F) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20230302_eAML_combined_mutated.png", width = 4, height = 4, dpi = 600, plot = p)


# genotyping rate of NPM1
length(intersect(eAML1.1.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML1.1"))
length(intersect(eAML1.2.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML1.2"))
length(intersect(eAML1.3.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML1.3"))
length(intersect(eAML1.4.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML1.4"))
length(intersect(eAML1.5.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML1.5"))
length(intersect(eAML1.6.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML1.6"))
length(intersect(eAML3.1.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML3.1"))
length(intersect(eAML3.2.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML3.2"))
length(intersect(eAML5.1.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML5.1"))
length(intersect(eAML5.2.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML5.2"))
length(intersect(eAML6.NPM1$bc, colnames(eAML.data))) / length(which(eAML.data$orig.ident == "eAML6"))
