# reanalysis of public data from Strobl et al., Br J Dermatology 2024

library(dplyr)
library(Seurat)

source("./scRNA_seq/R/00_eAML.colors.R")

for (sample in c("GvD_04", "GvD_05", "GvD_07", "GvD_08", "GvD_116")) {
  message(sample)
  seurat.data <- Read10X(paste0("./data/scRNA_seq/GSM7524599/", sample))
  HTO.data <- data.table::fread(paste0("./data/scRNA_seq/GSM7524599/", sample, "/HTO_demux.csv")) %>% as.data.frame()
  HTO.data$index <- paste0(HTO.data$index, "-1")
  rownames(HTO.data) <- HTO.data$index

  so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = sample, min.cells = 3, min.features = 200)
  so <- subset(so, cells = HTO.data$index)
  so$sample.id <- HTO.data[colnames(so), "hto_demux"]
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so <- RenameCells(so, add.cell.id = sample)
  if (sample == "GvD_04") {
    GvHD <- list(so)
  } else {
    GvHD <- c(GvHD, so)
  }
}

GvHD.combined <- merge(x = GvHD[[1]], y = list(GvHD[[2]], GvHD[[3]], GvHD[[4]], GvHD[[5]]), project = "GvHD")
GvHD.combined[["RNA"]] <- JoinLayers(GvHD.combined[["RNA"]])
GvHD.combined <- NormalizeData(GvHD.combined)
GvHD.combined <- FindVariableFeatures(GvHD.combined)
GvHD.combined <- ScaleData(GvHD.combined)
GvHD.combined <- RunPCA(GvHD.combined)
GvHD.combined <- RunUMAP(GvHD.combined, dims = 1:30)
GvHD.combined <- FindNeighbors(GvHD.combined)
GvHD.combined <- FindClusters(GvHD.combined, resolution = 0.2)

GvHD.combined <- SCTransform(GvHD.combined)

eAML.combined <- readRDS("./data/objects/20241021_eAML_combined.scvi.rds")
eAML.combined <- SCTransform(eAML.combined)

anchors <- FindTransferAnchors(
  reference = eAML.combined,
  query = GvHD.combined,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = 1:50
)
predictions <- TransferData(anchorset = anchors, refdata = eAML.combined$manual.cluster, dims = 1:30)
GvHD.combined <- AddMetaData(GvHD.combined, metadata = predictions)

saveRDS(file = "./data/scRNA_seq/objects/20250416_GvHD.rds", GvHD.combined)

GvHD.combined <- readRDS(file = "./data/scRNA_seq/objects/20250416_GvHD.rds")

GvHD.Tcells <- subset(GvHD.combined, predicted.id %in% c("T", "NK") & orig.ident != "GvD_116") # GvD_116 contains almost no T/NK cells

GvHD.Tcells <- NormalizeData(GvHD.Tcells)
GvHD.Tcells <- FindVariableFeatures(GvHD.Tcells)
GvHD.Tcells <- ScaleData(GvHD.Tcells)
GvHD.Tcells <- RunPCA(GvHD.Tcells)
GvHD.Tcells <- RunUMAP(GvHD.Tcells, dims = 1:30)
GvHD.Tcells <- FindNeighbors(GvHD.Tcells)
GvHD.Tcells <- FindClusters(GvHD.Tcells, resolution = 0.5)

markers <- FindAllMarkers(GvHD.Tcells, only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(GvHD.Tcells, features = top10$gene) + NoLegend()

GvHD.Tcells$manual.cluster <- "none"
GvHD.Tcells$manual.cluster[which(GvHD.Tcells$seurat_clusters %in% c("3", "4", "12"))] <- "CD8+ T cell"
GvHD.Tcells$manual.cluster[which(GvHD.Tcells$seurat_clusters %in% c("5"))] <- "prolif. CD8+ T cell"
GvHD.Tcells$manual.cluster[which(GvHD.Tcells$seurat_clusters %in% c("0", "1"))] <- "CD4+ T cell"
GvHD.Tcells$manual.cluster[which(GvHD.Tcells$seurat_clusters %in% c("6", "7"))] <- "NK"
GvHD.Tcells$manual.cluster[which(GvHD.Tcells$seurat_clusters %in% c("10"))] <- "prolif. NK cell"
GvHD.Tcells$manual.cluster[which(GvHD.Tcells$seurat_clusters %in% c("2", "8"))] <- "CD16+ NK"
saveRDS(GvHD.Tcells, file = "./data/scRNA_seq/objects/20250416_GvHD.Tcells.rds")
GvHD.Tcells <- readRDS(file = "./data/scRNA_seq/objects/20250416_GvHD.Tcells.rds")

GvHD.colors <- c(eAML.Tcell.colors, "none" = "grey90", "prolif. CD8+ T cell" = "firebrick", "prolif. NK cell" = "#164c07")

### overview of data
p <- DimPlot(subset(GvHD.Tcells, sample.id %in% c("HTO-skin", "HTO-blood")), group.by = "manual.cluster") +
  scale_color_manual(values = GvHD.colors) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20250505_GvHD_celltypes.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(subset(GvHD.Tcells, sample.id %in% c("HTO-skin", "HTO-blood")), group.by = "sample.id") +
  scale_color_manual(values = c("HTO-blood" = as.character(tissue.colors["PB"]), "HTO-skin" = "#92278F")) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20250505_GvHD_tissues.png", width = 4, height = 4, dpi = 600, plot = p)

df <- data.frame(
  cb = colnames(GvHD.Tcells),
  tissue = GvHD.Tcells$sample.id,
  sample = GvHD.Tcells$orig.ident,
  manual.cluster = GvHD.Tcells$manual.cluster
) %>%
  filter(tissue %in% c("HTO-skin", "HTO-blood")) %>%
  filter(manual.cluster != "none")

df.stats <- df %>%
  group_by(sample, manual.cluster, tissue) %>%
  summarize(cells = length(manual.cluster))

df.stats$cells.all <- apply(df.stats, MARGIN = 1, FUN = function(x) {
  length(which(GvHD.Tcells$orig.ident == x["sample"] &
    GvHD.Tcells$sample.id == x["tissue"] &
    GvHD.Tcells$manual.cluster != "none"))
})
df.stats$freq <- df.stats$cells / df.stats$cells.all
df.stats$sample.tissue <- paste0(df.stats$sample, ".", df.stats$tissue)
df.stats$manual.cluster <- factor(df.stats$manual.cluster, levels = c("CD16+ NK", "NK", "prolif. NK cell", "CD4+ T cell", "CD8+ T cell", "prolif. CD8+ T cell"))
ggplot(df.stats, aes(x = sample.tissue, y = 100 * freq)) +
  geom_col(aes(fill = manual.cluster), color = "black") +
  scale_fill_manual(values = GvHD.colors) +
  scale_x_discrete(labels = rep(c("blood", "skin"), 4)) +
  scale_y_continuous("% cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/GvHD/plots/20250506_overview_cell_freqs.svg", width = 1.5, height = 2.5)


### expression of exhaustion / memory score in blood versus GvHD
GvHD.Tcells <- AddModuleScore(GvHD.Tcells,
  features = list("Oliveira.tumor.score" = c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "TOX", "LAG3", "ENTPD1")),
  ctrl = 5, name = "Oliveira.tumor.score"
)
GvHD.Tcells <- AddModuleScore(GvHD.Tcells,
  features = list("Oliveira.memory.score" = c("TCF7", "IL7R", "SELL", "CCR7", "CD28")),
  ctrl = 5, name = "Oliveira.memory.score"
)

df <- data.frame(
  cb = colnames(GvHD.Tcells),
  tissue = GvHD.Tcells$sample.id,
  manual.cluster = GvHD.Tcells$manual.cluster,
  tumor.score = GvHD.Tcells$Oliveira.tumor.score1,
  memory.score = GvHD.Tcells$Oliveira.memory.score1
) %>%
  filter(tissue %in% c("HTO-skin", "HTO-blood"))

p <- ggplot(df[which(df$tissue == "HTO-blood" & df$manual.cluster == "CD4+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("Memory score", limits = c(-1, 2)) +
  scale_y_continuous("Exhaustion score", limits = c(-1, 2)) +
  NoLegend() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.text = element_text("Arial", size = 8, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20250505_exhaustion_memory_PB_GvHD_CD4.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "HTO-skin" & df$manual.cluster == "CD4+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("Memory score", limits = c(-1, 2)) +
  scale_y_continuous("Exhaustion score", limits = c(-1, 2)) +
  NoLegend() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.text = element_text("Arial", size = 8, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20250505_exhaustion_memory_skin_GvHD_CD4.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))


p <- ggplot(df[which(df$tissue == "HTO-blood" & df$manual.cluster == "prolif. CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("Memory score", limits = c(-1, 2)) +
  scale_y_continuous("Exhaustion score", limits = c(-1, 2)) +
  NoLegend() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.text = element_text("Arial", size = 8, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20250505_exhaustion_memory_PB_GvHD_CD8_proliferating.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "HTO-skin" & df$manual.cluster == "prolif. CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("Memory score", limits = c(-1, 2)) +
  scale_y_continuous("Exhaustion score", limits = c(-1, 2)) +
  NoLegend() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.text = element_text("Arial", size = 8, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20250505_exhaustion_memory_skin_GvHD_CD8_proliferating.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))


p <- ggplot(df[which(df$tissue == "HTO-blood" & df$manual.cluster == "CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("Memory score", limits = c(-1, 2)) +
  scale_y_continuous("Exhaustion score", limits = c(-1, 2)) +
  NoLegend() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.text = element_text("Arial", size = 8, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20250505_exhaustion_memory_PB_GvHD_CD8.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "HTO-skin" & df$manual.cluster == "CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_vline(xintercept = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("Memory score", limits = c(-1, 2)) +
  scale_y_continuous("Exhaustion score", limits = c(-1, 2)) +
  NoLegend() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.text = element_text("Arial", size = 8, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20250505_exhaustion_memory_skin_GvHD_CD8.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

df$exhausted <- ifelse(df$tumor.score > 0.5 & df$memory.score < 0.5, "exhausted", "non.exhausted")

df.stats <- df %>%
  group_by(manual.cluster, tissue) %>%
  summarize(exhausted.freq = length(which(exhausted == "exhausted")) /
    length(exhausted))

# proliferative fraction
df <- data.frame(
  cb = colnames(GvHD.Tcells),
  manual.cluster = GvHD.Tcells$manual.cluster,
  sample = GvHD.Tcells$orig.ident,
  tissue = GvHD.Tcells$sample.id
)

df.stats <- df %>%
  filter(tissue %in% c("HTO-skin", "HTO-blood")) %>%
  group_by(manual.cluster, sample, tissue) %>%
  summarize(cells = length(manual.cluster))

df.stats$cells.all <- sapply(df.stats$sample, FUN = function(x) {
  length(which(GvHD.Tcells$orig.ident == x))
})
df.stats$freq <- df.stats$cells / df.stats$cells.all



### differential gene expression between blood and GvHD

# CD8+ 
CD8.markers <- FindMarkers(subset(GvHD.Tcells, manual.cluster == "CD8+ T cell"), group.by = "sample.id", ident.1 = "HTO-blood", ident.2 = "HTO-skin")

genes.CD8 <- c(
  "CX3CR1", "GNLY", "PRF1", "NKG7", "GZMH", "GZMM", "GZMK", "KLRD1", "FCGR3A", "FCRL6", "ZAP70", "CXCL13", "IFNG", "ZFP36", "ZNF331", "CXCR4", "TNFRSF18", "ZFP36L1",
  "ICOS", "LAG3", "PDCD1", "CD44", "ITGB2", "KLRG1", "ITGAL", "ITGB7", "ITGA4"
)


genes.CD8 <- c(
  "CX3CR1", "CXCL13", "", "PRF1", "SELL", "ZFP36L1", "ZFP36", "ZFP36L2", "ZNF331", "CXCR4", "TNFRSF18", "ICOS",
  "LAG3", "PDCD1", "CD44", "ITGB2", "KLRG1", "ITGB7", "ITGA4"
)

CD8.markers$FDR <- -log10(CD8.markers$p_val_adj)
# CD8.markers = CD8.markers[which(abs(CD8.markers$avg_log2FC) < 2 & CD8.markers$FDR < 150),]
CD8.markers$gene <- rownames(CD8.markers)

ggplot() +
  ggrastr::rasterize(geom_point(data = CD8.markers, aes(x = avg_log2FC, y = FDR), color = "grey", size = 0.5), dpi = 600) +
  ggrastr::rasterize(geom_point(data = CD8.markers[which(CD8.markers$gene %in% genes.CD8), ], aes(x = avg_log2FC, y = FDR), color = "firebrick", size = 0.5), dpi = 600) +
  geom_label_repel(
    data = CD8.markers[which(CD8.markers$gene %in% genes.CD8), ], aes(x = avg_log2FC, y = FDR, label = gene),
    color = "black", label.size = 0, size = 2.5, max.overlaps = 20
  ) +
  scale_x_continuous("log2FC", limits = c(-13, 13)) +
  scale_y_continuous("-log10(FDR)", limits = c(0, 150)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/volcano/20250505_eAML_combined_CD8_PB_GvHD.svg", width = 4, height = 3)

# CD4+ 
CD4.markers <- FindMarkers(subset(GvHD.Tcells, manual.cluster == "CD4+ T cell"), group.by = "sample.id", ident.1 = "HTO-blood", ident.2 = "HTO-skin")

genes.CD4 <- c("ZFP36", "ZFP36L2", "CXCR4", "CTLA4", "ICOS", "CCL4", "CD44", "PDCD1", "ZNF331", "ITGA4", "SELL", "ITGB2", "TCF7", "ITGA6", "ITGB7", "ITGAL")

CD4.markers$FDR <- -log10(CD4.markers$p_val_adj)
# CD4.markers = CD4.markers[which(abs(CD4.markers$avg_log2FC) < 2 & CD4.markers$FDR < 150),]
CD4.markers$gene <- rownames(CD4.markers)

ggplot() +
  ggrastr::rasterize(geom_point(data = CD4.markers, aes(x = avg_log2FC, y = FDR), color = "grey", size = 0.5), dpi = 600) +
  ggrastr::rasterize(geom_point(data = CD4.markers[which(CD4.markers$gene %in% genes.CD4), ], aes(x = avg_log2FC, y = FDR), color = "firebrick", size = 0.5), dpi = 600) +
  geom_label_repel(
    data = CD4.markers[which(CD4.markers$gene %in% genes.CD4), ], aes(x = avg_log2FC, y = FDR, label = gene),
    color = "black", label.size = 0, size = 2.5, max.overlaps = 20
  ) +
  scale_x_continuous("log2FC", limits = c(-5, 5)) +
  scale_y_continuous("-log10(FDR)") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/volcano/20250505_eAML_combined_CD4_PB_GvHD.svg", width = 3, height = 3)
