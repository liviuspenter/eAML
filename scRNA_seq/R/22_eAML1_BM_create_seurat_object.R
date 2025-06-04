# create seurat object from BM data published in Maurer et al., Science Immunology 2025

library(dplyr)
library(ggplot2)
library(Seurat)

deconvolution.csv <- data.table::fread("./data/scRNA_seq/deconvolution/vireo_labels_allpools.csv") %>% as.data.frame()
deconvolution.csv$V1 <- gsub(deconvolution.csv$V1, pattern = "KMA1", replacement = "Pool97_23")
deconvolution.csv$V1 <- gsub(deconvolution.csv$V1, pattern = "KMA2", replacement = "Pool97_24")
deconvolution.csv$V1 <- gsub(deconvolution.csv$V1, pattern = "KMA3", replacement = "Pool97_25")
deconvolution.csv$V1 <- gsub(deconvolution.csv$V1, pattern = "KMA4", replacement = "Pool97_26")
deconvolution.csv$V1 <- gsub(deconvolution.csv$V1, pattern = "KMA5", replacement = "Pool97_27")
deconvolution.csv$V1 <- gsub(deconvolution.csv$V1, pattern = "KMA6", replacement = "Pool97_28")
deconvolution.csv$V1 <- paste0(deconvolution.csv$V1, "-1")
deconvolution.csv$sample <- substr(deconvolution.csv$V1, start = 1, stop = 9)
rownames(deconvolution.csv) <- deconvolution.csv$V1

for (s in c("Pool97_23", "Pool97_24", "Pool97_25", "Pool97_26", "Pool97_27", "Pool97_28")) {
  message(s)
  barcodes <- paste0(deconvolution.csv$barcode[which(deconvolution.csv$sample == s & deconvolution.csv$donor_id %in% c("donorLB", "tumorLB"))], "-1")

  seurat.data <- Read10X(data.dir = paste0("./data/scRNA_seq/", s, "/filtered_feature_bc_matrix/"))

  keep.cells <- which(colnames(seurat.data$`Gene Expression`) %in% barcodes)

  rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
  so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`[, keep.cells], project = s, min.cells = 3, min.features = 200)

  message(paste(s, " with ", ncol(so), " barcodes"))
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- subset(so, subset = nCount_RNA < 20000 & nFeature_RNA < 4000 & percent.mt < 20)
  message(paste(s, " with ", ncol(so), " barcodes after basic filtering"))

  so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
  so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")

  so <- RenameCells(so, add.cell.id = s)

  if (s == "Pool97_23") {
    AML44 <- list(so)
  } else {
    AML44 <- c(AML44, so)
  }
}

AML44 <- merge(AML44[[1]], list(AML44[[2]], AML44[[3]], AML44[[4]], AML44[[5]], AML44[[6]]))
AML44 <- JoinLayers(AML44)

AML44 <- NormalizeData(AML44)
AML44 <- FindVariableFeatures(AML44)
AML44 <- ScaleData(AML44)
AML44 <- RunPCA(AML44)
AML44 <- RunUMAP(AML44, dims = 1:30)
AML44 <- FindNeighbors(AML44)
AML44 <- FindClusters(AML44, resolution = 0.3)

# map to bone marrow reference - can be downloaded from Seurat
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
bm <- readRDS(file = "./data/scRNA_seq/objects/bm.reference.rds")
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./data/scRNA_seq/objects/reftmp.idx")

anchors <- FindTransferAnchors(
  reference = bm, query = AML44, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML44 <- MapQuery(
  anchorset = anchors, AML44, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)

AML44 <- AddMetaData(AML44, deconvolution.csv[, c("donor_id", "sample")])

# read TCR
source("./scRNA_seq/R/00_read_clonotypes.R")
source("./scRNA_seq/R/00_read_tcr_bcr.R")

metadata <- as.data.frame(readxl::read_excel("./data/scRNA_seq/20230207_TCR_eAML1_with_marrow.xlsx"))
rownames(metadata) <- metadata$Sample

tcr.data <- data.frame()
clonotypes <- data.frame()
for (s in metadata$Sample) {
  f <- read_tcr_bcr(paste0("./data/", metadata[s, "Pool"], "/all_contig_annotations.csv"),
    paste0("./data/", metadata[s, "Pool"], "/clonotypes.csv"),
    prefix = s, bcrtcr = "TCR"
  )
  tcr.data <- rbind(tcr.data, f)

  f <- data.table::fread(paste0("./data/", metadata[s, "Pool"], "/all_contig_annotations.csv")) %>% as.data.frame()
  f$barcode <- paste0(s, "_", f$barcode)
  f <- f %>% filter(barcode %in% colnames(AML44) & is_cell == TRUE)
  if (nrow(f) == 0) {
    next()
  }
  clonotypes <- rbind(clonotypes, f, fill = T)
}
AML44 <- AddMetaData(AML44, metadata = tcr.data)
AML44$TCR_clonotype_id <- factor(AML44$TCR_clonotype_id, levels = gtools::mixedsort(unique(AML44$TCR_clonotype_id)))

saveRDS(AML44, file = "./data/scRNA_seq/objects/AML44_BM.rds")
