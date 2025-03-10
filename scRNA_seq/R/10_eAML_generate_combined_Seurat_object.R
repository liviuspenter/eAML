# generate combined Seurat object of BM, PB, skin and eAML samples
# integration with scvi tools

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(Seurat)

source("./scRNA_seq/analysis/00_eAML.colors.R")

# read metadata
metadata <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 3))
rownames(metadata) <- metadata$sample

eAML.genotype <- data.frame()
for (s in metadata$sample) {
  f <- as.data.frame(data.table::fread(paste0(metadata[s, "library"], "/souporcell/clusters.tsv")))
  f$barcode <- paste0(s, "_", f$barcode)
  eAML.genotype <- rbind(eAML.genotype, f)
}
rownames(eAML.genotype) <- eAML.genotype$barcode

metadata <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 1))
rownames(metadata) <- metadata$sample

# create individual Seurat objects
so.list <- list()
for (s in metadata$sample) {
  message(s)
  seurat.data <- Read10X(metadata[s, "library"])
  if (length(seurat.data) == 2) {
    eAML.data <- CreateSeuratObject(seurat.data$`Gene Expression`, project = s, min.cells = 3, min.features = 200)
  } else {
    eAML.data <- CreateSeuratObject(seurat.data, project = s, min.cells = 3, min.features = 200)
  }

  eAML.data[["percent.mt"]] <- PercentageFeatureSet(eAML.data, pattern = "^MT-")
  eAML.data <- RenameCells(eAML.data, add.cell.id = s)
  eAML.data$patient <- metadata[s, "patient"]
  eAML.data$tissue <- metadata[s, "tissue"]

  if (s == metadata$sample[1]) {
    so.list <- eAML.data
  } else {
    so.list <- c(so.list, eAML.data)
  }
}
eAML.data <- merge(x = so.list[[1]], y = so.list[2:length(so.list)], project = "eAML.data")

# exclude doublets samples where such data is available
eAML.data <- subset(eAML.data, cells = c(
  colnames(eAML.data)[which(!eAML.data$orig.ident %in%
    unique(stringr::str_split_fixed(eAML.genotype$barcode, pattern = "_", n = 2)[, 1]))],
  eAML.genotype$barcode[which(eAML.genotype$status == "singlet")]
))

# read TCR
source("./scRNA_seq/R/00_read_clonotypes.R")
source("./scRNA_seq/R/00_read_tcr_bcr.R")

metadata <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 2))
rownames(metadata) <- metadata$sample

tcr.data <- data.frame()
clonotypes <- data.frame()
for (s in metadata$sample) {
  f <- read_tcr_bcr(paste0(metadata[s, "library"], "/all_contig_annotations.csv"),
    paste0(metadata[s, "library"], "clonotypes.csv"),
    prefix = s, bcrtcr = "TCR"
  )
  tcr.data <- rbind(tcr.data, f)

  f <- data.table::fread(paste0(metadata[s, "library"], "/all_contig_annotations.csv"))
  f$barcode <- paste0(s, "_", f$barcode)
  f <- f %>% filter(barcode %in% colnames(eAML.data) & is_cell == TRUE)
  clonotypes <- rbind(clonotypes, f)
}
eAML.data <- AddMetaData(eAML.data, metadata = tcr.data)
eAML.data$TCR_clonotype_id <- factor(eAML.data$TCR_clonotype_id, levels = gtools::mixedsort(unique(eAML.data$TCR_clonotype_id)))

eAML.data <- NormalizeData(eAML.data)
eAML.data <- FindVariableFeatures(eAML.data)
eAML.data <- ScaleData(eAML.data)
eAML.data <- RunPCA(eAML.data)
eAML.data <- RunUMAP(eAML.data, dims = 1:30)
eAML.data <- FindNeighbors(eAML.data)
eAML.data <- FindClusters(eAML.data, resolution = 0.4)
DimPlot(eAML.data, group.by = "tissue", cols = tissue.colors)
saveRDS(file = "./data/scRNA_seq/objects/20240802_eAML_combined.rds", eAML.data)

eAML.data <- readRDS(file = "./data/scRNA_seq/objects/20240802_eAML_combined.rds")

# integration with scvi

library(reticulate)
library(anndata)
library(sceasy)

# need to adjust to local conda environment
use_condaenv("/Users/liviuspenter/miniconda3/envs/scvi-env/")

sc <- import("scanpy")
scvi <- import("scvi")

eAML.data[["RNA"]] <- as(object = eAML.data[["RNA"]], Class = "Assay")
adata <- convertFormat(eAML.data, from = "seurat", to = "anndata", main_layer = "counts", drop_single_values = FALSE)
scvi$model$SCVI$setup_anndata(adata, batch_key = "orig.ident")
model <- scvi$model$SCVI(adata)
# this step takes a long time
model$train()
latent <- model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) <- colnames(eAML.data)
eAML.data[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(eAML.data))
eAML.data <- FindNeighbors(eAML.data, reduction = "scvi")
eAML.data <- FindClusters(eAML.data, resolution = 0.1)
eAML.data <- RunUMAP(eAML.data, dims = 1:10, reduction = "scvi")
saveRDS(file = "./data/scRNA_seq/objects/20241021_eAML_combined.scvi.rds", eAML.data)
eAML.data <- readRDS(file = "./data/scRNA_seq/objects/20241021_eAML_combined.scvi.rds")

# annotate clusters
eAML.data.markers <- FindAllMarkers(eAML.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(file = "./data/scRNA_seq/objects/20241021_eAML_combined_markers.scvi.rds", eAML.data.markers)
View(eAML.data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC))

markers <- as.data.frame(eAML.data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  select(gene))
markers$cluster <- factor(markers$cluster, levels = names(eAML.combined.colors))
markers <- with(markers, markers[order(cluster), ])

# visualize marker genes
boo <- subset(eAML.data, downsample = 400)
boo <- ScaleData(boo, features = unique(markers$gene))
boo$manual.cluster <- factor(boo$manual.cluster, levels = names(eAML.combined.colors))
p <- DoHeatmap(boo,
  group.by = "manual.cluster", features = unique(markers$gene), group.bar.height = 0.01,
  group.colors = eAML.combined.colors, label = F
) + NoLegend() +
  scale_fill_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")), na.value = "black") +
  theme(axis.text = element_text("Arial", size = 6, color = "black", face = "italic"))
ggsave("./scRNA_seq/figures/combined/heatmap/20241016_eAML_combined_marker_genes.png", width = 5, height = 10, plot = p, dpi = 600)

eAML.data$manual.cluster <- "none"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(0))] <- "T"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(3))] <- "NK"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(1, 6, 7, 15, 17))] <- "myeloid"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(5, 16))] <- "B cell"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(2))] <- "macrophage"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(4))] <- "fibroblast"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(10))] <- "erythroid"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(8))] <- "keratinocyte"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(14))] <- "platelet"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(9, 11))] <- "endothelium"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(13))] <- "melanocyte"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(12))] <- "smooth muscle"

# integrate chimerism data
AML1007.genotype <- as.data.frame(read.table("./data/scRNA_seq/souporcell/20210707_AML1007_chimerism.csv"))
AML1007.genotype <- AML1007.genotype[which(grepl("AML1007.1", AML1007.genotype$barcode)), ]
AML1007.genotype$barcode <- sub(AML1007.genotype$barcode, pattern = "AML1007.1", replacement = "1007")
AML1010.genotype <- read.table("./data/scRNA_seq/souporcell/20210707_AML1010_chimerism.csv")
AML1010.genotype <- AML1010.genotype[which(grepl("AML1010.1", AML1010.genotype$barcode)), ]
AML1010.genotype$barcode <- sub(AML1010.genotype$barcode, pattern = "AML1010.1", replacement = "1010")
AML1012.genotype <- read.table("./data/scRNA_seq/souporcell/20210707_AML1012_chimerism.csv")
AML1012.genotype <- AML1012.genotype[which(grepl("AML1012.1", AML1012.genotype$barcode)), ]
AML1012.genotype$barcode <- sub(AML1012.genotype$barcode, pattern = "AML1012.1", replacement = "1012")
AML1016.genotype <- read.table("./data/scRNA_seq/souporcell/20210707_AML1016_chimerism.csv")
AML1016.genotype <- AML1016.genotype[which(grepl("AML1016.1", AML1016.genotype$barcode)), ]
AML1016.genotype$barcode <- sub(AML1016.genotype$barcode, pattern = "AML1016.1", replacement = "1016")
AML1019.genotype <- read.table("./data/scRNA_seq/souporcell/20210707_AML1019_chimerism.csv")
AML1019.genotype <- AML1019.genotype[which(grepl("AML1019.1", AML1019.genotype$barcode)), ]
AML1019.genotype$barcode <- sub(AML1019.genotype$barcode, pattern = "AML1019.1", replacement = "1019")
AML1022.genotype <- read.table("./data/scRNA_seq/souporcell/20210707_AML1022_chimerism.csv")
AML1022.genotype <- AML1022.genotype[which(grepl("AML1022.1", AML1022.genotype$barcode)), ]
AML1022.genotype$barcode <- sub(AML1022.genotype$barcode, pattern = "AML1022.1", replacement = "1022")
AML1026.genotype <- read.table("./data/scRNA_seq/souporcell/20210707_AML1026_chimerism.csv")
AML1026.genotype <- AML1026.genotype[which(grepl("AML1026.1", AML1026.genotype$barcode)), ]
AML1026.genotype$barcode <- sub(AML1026.genotype$barcode, pattern = "AML1026.1", replacement = "1026")
eAML1.genotype <- read.table("./data/scRNA_seq/souporcell/20230119_eAML1_genotype.csv")
eAML2.genotype <- read.table("./data/scRNA_seq/souporcell/20230119_eAML2_genotype.csv")
eAML4.genotype <- read.table("./data/scRNA_seq/souporcell/20230119_eAML4_genotype.all.csv")
eAML5.genotype <- read.table("./data/scRNA_seq/souporcell/20230119_eAML5_genotype.csv")
eAML6.genotype <- read.table("./data/scRNA_seq/souporcell/20230119_eAML6_genotype.csv")
eAML7.genotype <- read.table("./data/scRNA_seq/souporcell/20230119_eAML7_genotype.csv")

AML.genotype <- as.data.frame(dplyr::bind_rows(
  AML1007.genotype, AML1010.genotype, AML1012.genotype, AML1016.genotype,
  AML1019.genotype, AML1022.genotype, AML1026.genotype,
  eAML1.genotype, eAML2.genotype, eAML4.genotype, eAML5.genotype, eAML6.genotype, eAML7.genotype
))
rownames(AML.genotype) <- AML.genotype$barcode
eAML.data$chimerism <- AML.genotype[colnames(eAML.data), "chimerism"]

# plot UMAPs
p <- DimPlot(eAML.data, shuffle = T, group.by = "chimerism", cols = donor.recipient.colors, na.value = "grey90", raster = F) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241014_eAML_combined_donor_recipient.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(eAML.data, shuffle = T, group.by = "manual.cluster", cols = eAML.combined.colors, na.value = "grey90", raster = F) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241014_eAML_combined_manual_clusters.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(eAML.data, shuffle = T, group.by = "tissue", cols = tissue.colors, na.value = "grey90", raster = F) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241014_eAML_combined_tissue.png", width = 4, height = 4, dpi = 600, plot = p)

# visualize celltype fractions
Idents(eAML.data) <- "manual.cluster"
prop.data <- as.data.frame(prop.table(table(Idents(eAML.data), eAML.data$orig.ident), margin = 2))
prop.data$Var2 <- factor(prop.data$Var2, levels = c(BM.samples, PB.samples, eAML.samples2, skin.samples))
prop.data$Var1 <- factor(prop.data$Var1, levels = names(eAML.combined.colors))
ggplot(prop.data, aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col(color = "black") +
  scale_fill_manual(values = eAML.combined.colors) +
  scale_y_continuous("% cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241015_overview_celltypes.svg", width = 4, height = 2.5)

# statistics of celltype frequencies
col_fun <- circlize::colorRamp2(breaks = seq(-2, 2, 4 / 8), colors = BuenColors::jdb_palette(name = "brewer_celsius"))
prop.data <- as.data.frame.matrix(prop.table(table(Idents(eAML.data), eAML.data$orig.ident), margin = 2))
prop.data <- prop.data[, c(BM.samples, PB.samples, eAML.samples, skin.samples)]
prop.data$manual.cluster <- rownames(prop.data)
prop.data <- as.data.frame(prop.data %>% tidyr::pivot_longer(cols = colnames(prop.data)[1:ncol(prop.data) - 1], names_to = "sample"))
prop.data$sample <- factor(prop.data$sample, levels = c(BM.samples, PB.samples, eAML.samples, skin.samples))
prop.data$tissue <- "BM"
prop.data$tissue[which(prop.data$sample %in% PB.samples)] <- "PB"
prop.data$tissue[which(prop.data$sample %in% eAML.samples)] <- "eAML"
prop.data$tissue[which(prop.data$sample %in% skin.samples)] <- "skin"
prop.data$tissue <- factor(prop.data$tissue, levels = c("BM", "PB", "eAML", "skin"))

ggplot(prop.data[which(prop.data$manual.cluster == "myeloid"), ], aes(x = tissue, y = 100 * as.numeric(value))) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = names(tissue.colors)) +
  scale_y_continuous("% progenitor") +
  scale_color_manual(values = tissue.colors) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 1) +
  geom_signif(comparisons = list(c("PB", "eAML"), c("BM", "PB")), step_increase = 0.1, textsize = 3, test = "t.test") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241015_eAML_combined_myeloid.svg", width = 1.5, height = 2.2)

ggplot(prop.data[which(prop.data$manual.cluster == "T"), ], aes(x = tissue, y = 100 * as.numeric(value))) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_point(aes(color = tissue), size = 0.5) +
  scale_x_discrete(limits = names(tissue.colors)) +
  scale_y_continuous("% T cells") +
  scale_color_manual(values = tissue.colors) +
  geom_signif(comparisons = list(c("BM", "PB"), c("PB", "eAML")), step_increase = 0.1, textsize = 3, test = "t.test") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241015_eAML_combined_T.svg", width = 1.3, height = 2)

ggplot(prop.data[which(prop.data$manual.cluster == "NK"), ], aes(x = tissue, y = 100 * as.numeric(value))) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_point(aes(color = tissue), size = 0.5) +
  scale_x_discrete(limits = names(tissue.colors)) +
  scale_y_continuous("% NK cells") +
  scale_color_manual(values = tissue.colors) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("PB", "eAML")), step_increase = 0.1, textsize = 3, test = "t.test") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241015_eAML_combined_NK.svg", width = 1.3, height = 2)


df <- data.frame(
  sample = eAML.data$orig.ident,
  tissue = eAML.data$tissue,
  manual.cluster = eAML.data$manual.cluster
)

df <- df %>%
  filter(manual.cluster %in% c("myeloid", "T", "NK", "B cell", "macrophage", "erythroid", "platelet")) %>%
  group_by(sample, tissue, manual.cluster) %>%
  tally()
df$total <- apply(df, MARGIN = 1, FUN = function(x) {
  sum(df[(which(df$sample == x["sample"])), "n"])
})
df$freq <- df$n / df$total

ggplot(df[which(df$manual.cluster == "T"), ], aes(x = tissue, y = 100 * as.numeric(freq))) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = names(tissue.colors)) +
  scale_y_continuous("% T cells (immune cells)") +
  scale_color_manual(values = tissue.colors) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("skin", "eAML")), step_increase = 0.1, textsize = 3, test = "t.test") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241022_eAML_combined_T.svg", width = 1.5, height = 2.2)


ggplot(df[which(df$manual.cluster == "NK"), ], aes(x = tissue, y = 100 * as.numeric(freq))) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = names(tissue.colors)) +
  scale_y_continuous("% NK cells (immune cells)") +
  scale_color_manual(values = tissue.colors) +
  geom_signif(comparisons = list(c("PB", "eAML")), step_increase = 0.1, textsize = 3, test = "t.test") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241022_eAML_combined_NK.svg", width = 1.5, height = 2.2)


stats <- df %>%
  group_by(sample, tissue) %>%
  summarize(ratio = n[which(manual.cluster == "T")] / n[which(manual.cluster == "NK")])

ggplot(stats, aes(x = tissue, y = as.numeric(ratio))) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = names(tissue.colors)) +
  scale_y_log10("ratio T / NK cells") +
  scale_color_manual(values = tissue.colors) +
  geom_signif(comparisons = list(c("BM", "eAML")), step_increase = 0.1, textsize = 3, test = "t.test") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241022_eAML_combined_ratio_T_NK.svg", width = 1.5, height = 2.2)

Idents(eAML.data) <- "manual.cluster"
col_fun <- circlize::colorRamp2(breaks = seq(-2, 2, 4 / 8), colors = BuenColors::jdb_palette(name = "brewer_celsius"))
prop.data <- prop.table(table(Idents(eAML.data), eAML.data$orig.ident), margin = 2)
prop.data <- as.data.frame.matrix(t(scale(t(prop.data))))

ha <- HeatmapAnnotation(
  tissue = factor(c(rep("BM", 7), rep("eAML", 16), rep("PB", 6), rep("skin", 6)), levels = c("BM", "PB", "eAML", "skin")),
  col = list("tissue" = c("BM" = "firebrick", "PB" = "red", "eAML" = "orange", "skin" = "blue")), border = T,
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10), simple_anno_size = unit(5, "pt")
)
svglite::svglite("./scRNA_seq/figures/combined/heatmap/20241016_eAML_combined_celltypes.svg", width = 4, height = 2)
Heatmap(
  as.matrix(prop.data[c(
    "myeloid", "macrophage", "T", "NK", "B cell", "fibroblast", "keratinocyte", "endothelium", "melanocyte",
    "smooth muscle", "erythroid", "platelet"
  ), ]),
  col = col_fun,
  cluster_columns = F, cluster_rows = F, border = T, column_title_gp = gpar(fontsize = 10),
  column_split = factor(c(rep("BM", 7), rep("eAML", 16), rep("PB", 6), rep("skin", 6)), levels = c("BM", "PB", "eAML", "skin")),
  top_annotation = ha, show_column_names = F, row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  row_labels = c(
    "progenitor", "mono/macropahge", "T", "NK", "B cell", "fibroblast", "keratinocyte", "endothelium", "melanocyte",
    "smooth muscle", "erythroid", "megakaryocytic"
  )
)
dev.off()

# tissue-specific chimerism
df <- data.frame(
  barcode = colnames(eAML.data),
  manual.cluster = eAML.data$manual.cluster,
  tissue = eAML.data$tissue,
  chimerism = eAML.data$chimerism,
  sample = eAML.data$orig.ident
)

df <- as.data.frame(df %>% filter(tissue != "skin") %>%
  filter(!is.na(chimerism)) %>%
  group_by(sample, manual.cluster) %>%
  filter(length(sample) > 5) %>%
  summarize(
    donor = length(which(chimerism == "donor")),
    recipient = length(which(chimerism == "recipient")),
    tissue = unique(tissue),
    sample = unique(sample)
  ))

df <- df[-which((df$donor + df$recipient) < 6), ]

df$donor.chimerism <- df$donor / (df$donor + df$recipient)

df <- df[-which(df$manual.cluster %in% c("erythroid", "platelet")), ]
df$manual.cluster.tissue <- paste0(df$manual.cluster, ".", df$tissue)
df <- df[-which(df$manual.cluster.tissue %in% c("keratinocyte.BM", "fibroblast.BM")), ]
df$manual.cluster.tissue <- factor(df$manual.cluster.tissue, levels = c(
  "endothelium.eAML", "fibroblast.eAML", "keratinocyte.eAML", "melanocyte.eAML",
  "smooth muscle.eAML",
  "B cell.BM", "B cell.PB", "B cell.eAML",
  "macrophage.BM", "macrophage.PB", "macrophage.eAML",
  "myeloid.BM", "myeloid.PB", "myeloid.eAML",
  "T.BM", "T.PB", "T.eAML",
  "NK.BM", "NK.PB", "NK.eAML"
))

ggplot(df, aes(x = manual.cluster.tissue, y = 100 * donor.chimerism, group = tissue, color = manual.cluster)) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.8, color = "black") +
  geom_point(size = 0.5) +
  scale_color_manual(values = eAML.combined.colors) +
  scale_y_continuous("% donor chimerism") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241016_donor_chimerism.svg", width = 3.5, height = 2.5)


ggplot(df, aes(x = manual.cluster.tissue, y = 100 * donor.chimerism, group = tissue, color = manual.cluster)) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.8, color = "black") +
  geom_point(size = 1) +
  scale_color_manual(values = eAML.combined.colors) +
  scale_y_continuous("% donor chimerism") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20250227_donor_chimerism.svg", width = 4.5, height = 2.5)

# more specialized analyses into tissue-specific chimerism or MRD
MRD.df <- data.frame(
  cb = colnames(eAML.data),
  tissue = eAML.data$tissue,
  chimerism = eAML.data$chimerism,
  sample = eAML.data$orig.ident,
  manual.cluster = eAML.data$manual.cluster
)

stats <- MRD.df[which(MRD.df$chimerism %in% c("donor", "recipient") & MRD.df$tissue == "PB"), ] %>%
  group_by(sample, manual.cluster) %>%
  summarize(
    donor = length(which(chimerism == "donor")),
    recipient = length(which(chimerism == "recipient"))
  )
stats$total <- sapply(stats$sample, FUN = function(x) {
  length(which(eAML.data$orig.ident == x))
})
stats$MRD.celltype <- stats$recipient / (stats$recipient + stats$donor)
stats$MRD.total <- stats$recipient / stats$total


stats$sample <- factor(stats$sample, levels = c("eAML2PB", "eAML4PB.2", "eAML4PB.1", "eAML4PB.3", "eAML7PB"))

# eAML4PB.2 - 7/22/22 - clinical T cell chimerism 34%
# eAML4PB.1 - 8/29/22 - clinical T cell chimerism 46%
# eAML4PB.4 - 1/27/23 - clinical T cell chimerism 100%

ggplot(stats[which(stats$manual.cluster == "T" & stats$sample %in% c("eAML4PB.1", "eAML4PB.2", "eAML4PB.3")), ], aes(x = sample, y = 100 * MRD.celltype)) +
  geom_col() +
  theme_classic() +
  scale_x_discrete(labels = c("pre-ipi", "pre-ipi", "post-ipi")) +
  scale_y_continuous("% recipient T cells") +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/eAML4/20241213_eAML4_Tcell_chimerism.svg", width = 1.3, height = 2)

length(which(eAML.data$orig.ident == "eAML4.1" & eAML.data$chimerism == "recipient" & eAML.data$manual.cluster == "T")) /
  length(which(eAML.data$orig.ident == "eAML4.1" & eAML.data$manual.cluster == "T"))

length(which(eAML.data$orig.ident == "eAML4.2" & eAML.data$chimerism == "recipient" & eAML.data$manual.cluster == "T")) /
  length(which(eAML.data$orig.ident == "eAML4.2" & eAML.data$manual.cluster == "T"))


ggplot(stats[which(stats$manual.cluster == "myeloid" & stats$sample %in% c("eAML2PB", "eAML7PB")), ], aes(x = sample, y = 100 * MRD.celltype)) +
  geom_col(fill = "#d62728", color = "black") +
  theme_classic() +
  scale_x_discrete(labels = c("eAML2", "eAML7")) +
  scale_y_continuous("% recipient myeloid cells") +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241213_eAML_MRD_myeloid.svg", width = 1.1, height = 2)

ggplot(stats[which(stats$manual.cluster == "myeloid" & stats$sample %in% c("eAML2PB", "eAML7PB")), ], aes(x = sample, y = 100 * MRD.total)) +
  geom_col(fill = "grey", color = "black") +
  theme_classic() +
  scale_x_discrete(labels = c("eAML2", "eAML7")) +
  scale_y_continuous("% recipient cells") +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241213_eAML_MRD_all.svg", width = 1.1, height = 2)
