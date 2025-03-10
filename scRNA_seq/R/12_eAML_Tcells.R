# focused analyis of T cells

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(Seurat)

source("./analysis/00_eAML.colors.R")

eAML.data <- readRDS("./data/objects/20241021_eAML_combined.scvi.rds")

eAML.data.Tcells <- subset(eAML.data, manual.cluster %in% c("T", "NK"))

library(reticulate)
library(anndata)
library(sceasy)

# adapt to local environment
use_condaenv("/Users/liviuspenter/miniconda3/envs/scvi-env/")

sc <- import("scanpy")
scvi <- import("scvi")

eAML.data.Tcells[["RNA"]] <- as(object = eAML.data.Tcells[["RNA"]], Class = "Assay")
adata <- convertFormat(eAML.data.Tcells, from = "seurat", to = "anndata", main_layer = "counts", drop_single_values = FALSE)

scvi$model$SCVI$setup_anndata(adata, batch_key = "orig.ident")
model <- scvi$model$SCVI(adata)
model$train()
latent <- model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) <- colnames(eAML.data.Tcells)
eAML.data.Tcells[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(eAML.data.Tcells))
eAML.data.Tcells <- FindNeighbors(eAML.data.Tcells, reduction = "scvi")
eAML.data.Tcells <- FindClusters(eAML.data.Tcells, resolution = 0.2)
eAML.data.Tcells <- RunUMAP(eAML.data.Tcells, dims = 1:10, reduction = "scvi", n.components = 2)
saveRDS(eAML.data.Tcells, file = "./data/scRNA_seq/objects/20241017_eAML_combined_Tcells.rds")

markers <- FindAllMarkers(eAML.data.Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC))

eAML.data.Tcells$manual.cluster <- "none"
eAML.data.Tcells$manual.cluster[which(eAML.data.Tcells$seurat_clusters %in% c(0, 7))] <- "CD4+ T cell"
eAML.data.Tcells$manual.cluster[which(eAML.data.Tcells$seurat_clusters %in% c(2, 9, 10, 11))] <- "CD8+ T cell"
eAML.data.Tcells$manual.cluster[which(eAML.data.Tcells$seurat_clusters %in% c(5))] <- "NK"
eAML.data.Tcells$manual.cluster[which(eAML.data.Tcells$seurat_clusters %in% c(1))] <- "CD16+ NK"
eAML.data.Tcells$manual.cluster[which(eAML.data.Tcells$seurat_clusters %in% c(6))] <- "Treg"
eAML.data.Tcells$manual.cluster[which(eAML.data.Tcells$seurat_clusters %in% c(4))] <- "naive_CM"
eAML.data.Tcells$manual.cluster[which(eAML.data.Tcells$seurat_clusters %in% c(3, 8, 12))] <- "doublet"
saveRDS(file = "./data/scRNA_seq/objects/20241017_eAML_combined_Tcells.rds", eAML.data.Tcells)
eAML.data.Tcells <- readRDS("./data/scRNA_seq/objects/20241017_eAML_combined_Tcells.rds")

clusters.markers <- FindAllMarkers(subset(eAML.data.Tcells, manual.cluster != "doublet"), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- as.data.frame(clusters.markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC) %>%
  dplyr::select(gene))
markers$cluster <- factor(markers$cluster, levels = names(eAML.combined.colors))
markers <- with(markers, markers[order(cluster), ])

boo <- subset(eAML.data.Tcells, downsample = 400)
boo <- subset(boo, manual.cluster != "doublet")
boo <- ScaleData(boo, features = unique(markers$gene))
boo$manual.cluster <- factor(boo$manual.cluster, levels = names(eAML.Tcell.colors))

p <- DoHeatmap(boo,
  group.by = "manual.cluster", features = unique(markers$gene), group.bar.height = 0.01,
  group.colors = eAML.Tcell.colors, label = F
) + NoLegend() +
  scale_fill_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")), na.value = "black") +
  theme(axis.text = element_text("Arial", size = 6, color = "black", face = "italic"))
ggsave("./scRNA_seq/figures/combined/heatmap/20241216_eAML_Tcells_combined_marker_genes.png", width = 4, height = 6.5, plot = p, dpi = 600)


p <- DimPlot(eAML.data.Tcells, group.by = "manual.cluster", cols = eAML.Tcell.colors) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241023_eAML_combined_Tcells_clusters.png", width = 4, height = 4, dpi = 600, plot = p)

# statistics of T/NK cell subsets
boo <- subset(eAML.data.Tcells, manual.cluster != "doublet")
Idents(eAML.data.Tcells) <- "manual.cluster"
prop.data <- as.data.frame(prop.table(table(Idents(boo), boo$orig.ident), margin = 2))
prop.data$Var2 <- factor(prop.data$Var2, levels = c(BM.samples, PB.samples, eAML.samples2, skin.samples))
prop.data$Var1 <- factor(prop.data$Var1, levels = names(eAML.Tcell.colors))

ggplot(prop.data, aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col(color = "black") +
  scale_fill_manual(values = eAML.Tcell.colors) +
  scale_y_continuous("% cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_TNK_cells_overview_celltypes.svg", width = 4, height = 2.5)

prop.data <- prop.table(table(Idents(boo), boo$orig.ident), margin = 2)
prop.data <- as.data.frame.matrix(t(prop.data))
prop.data <- prop.data[c(BM.samples, PB.samples, eAML.samples, skin.samples), ]
prop.data$tissue <- c(
  rep("BM", length(BM.samples)),
  rep("PB", length(PB.samples)),
  rep("eAML", length(eAML.samples)),
  rep("skin", length(skin.samples))
)

ggplot(prop.data, aes(x = tissue, y = 100 * Treg)) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = c("BM", "PB", "eAML", "skin")) +
  scale_y_continuous("% Treg") +
  scale_color_manual(values = tissue.colors) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("PB", "eAML"), c("BM", "skin"), c("PB", "skin")), textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_eAML_Tregs.svg", width = 1.1, height = 1.5)

ggplot(prop.data, aes(x = tissue, y = 100 * naive_CM)) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = c("BM", "PB", "eAML", "skin")) +
  scale_y_continuous("% naive / CM") +
  scale_color_manual(values = tissue.colors) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("PB", "eAML"), c("BM", "skin"), c("PB", "skin")), textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_eAML_naive_CM.svg", width = 1.1, height = 1.5)

ggplot(prop.data, aes(x = tissue, y = 100 * `CD16+ NK`)) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = c("BM", "PB", "eAML", "skin")) +
  scale_y_continuous("% CD16+ NK cells") +
  scale_color_manual(values = tissue.colors) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("BM", "skin")), textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_eAML_CD16_NK.svg", width = 1.1, height = 1.5)

ggplot(prop.data, aes(x = tissue, y = 100 * `CD4+ T cell`)) +
  geom_point(aes(color = tissue), size = 0.5) +
  scale_x_discrete(limits = c("BM", "PB", "eAML", "skin")) +
  scale_y_continuous("% CD4+ T cell") +
  scale_color_manual(values = tissue.colors) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_signif(comparisons = list(c("eAML", "skin")), textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_eAML_CD4_Tcell.svg", width = 1.1, height = 1.5)

ggplot(prop.data, aes(x = tissue, y = 100 * `CD8+ T cell`)) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_discrete(limits = c("BM", "PB", "eAML", "skin")) +
  scale_y_continuous("% CD8+ T cell") +
  scale_color_manual(values = tissue.colors) +
  stat_summary(geom = "crossbar", fun = median, fun.max = median, fun.min = median, size = 0.5, width = 0.5) +
  geom_signif(comparisons = list(c("eAML", "skin")), textsize = 3, step_increase = 0.1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_eAML_CD8_Tcell.svg", width = 1.1, height = 1.5)

ggplot(prop.data, aes(x = 100 * Treg, y = 100 * naive_CM)) +
  geom_point(aes(color = tissue), size = 1) +
  scale_x_continuous("% Tregs") +
  scale_y_continuous("% naive / CM") +
  scale_color_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_eAML_naive_CM_vs_Tregs.svg", width = 1.3, height = 1.3)

ggplot(prop.data, aes(x = 100 * `CD4+ T cell`, y = 100 * `CD8+ T cell`)) +
  geom_point(aes(color = tissue), size = 0.5) +
  scale_x_continuous("% Tregs") +
  scale_y_continuous("% naive / CM") +
  scale_color_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20241023_eAML_naive_CM_vs_Tregs.svg", width = 1.3, height = 1.3)


eAML.data.Tcells.CD8 <- subset(eAML.data.Tcells, cells = colnames(eAML.data.Tcells)[which(eAML.data.Tcells$manual.cluster == "CD8+ T cell")])
DefaultAssay(eAML.data.Tcells.CD8) <- "RNA"
eAML.data.Tcells.CD8 <- NormalizeData(eAML.data.Tcells.CD8)
eAML.data.Tcells.CD8 <- ScaleData(eAML.data.Tcells.CD8)
eAML.data.Tcells.CD8 <- FindVariableFeatures(eAML.data.Tcells.CD8)
eAML.data.Tcells.CD8 <- RunPCA(eAML.data.Tcells.CD8)
eAML.data.Tcells.CD8 <- RunUMAP(eAML.data.Tcells.CD8, dims = 1:10)

CD8.markers <- FindMarkers(eAML.data.Tcells.CD8, group.by = "tissue", ident.1 = "BM", ident.2 = "eAML")

genes.CD8 <- c(
  "CX3CR1", "GNLY", "PRF1", "NKG7", "GZMH", "GZMM", "GZMK", "KLRD1", "FCGR3A", "FCRL6", "ZAP70", "CXCL13", "IFNG", "ZFP36", "ZNF331", "CXCR4", "TNFRSF18", "ZFP36L1",
  "ICOS", "LAG3", "PDCD1", "CD44", "ITGB2", "KLRG1", "ITGAL", "ITGB7", "ITGA4"
)

CD8.markers$FDR <- -log10(CD8.markers$p_val_adj)
CD8.markers <- CD8.markers[which(abs(CD8.markers$avg_log2FC) < 2 & CD8.markers$FDR < 200), ]
CD8.markers$gene <- rownames(CD8.markers)

ggplot() +
  ggrastr::rasterize(geom_point(data = CD8.markers, aes(x = avg_log2FC, y = FDR), color = "grey", size = 0.5), dpi = 600) +
  ggrastr::rasterize(geom_point(data = CD8.markers[which(CD8.markers$gene %in% genes.CD8), ], aes(x = avg_log2FC, y = FDR), color = "firebrick", size = 0.5), dpi = 600) +
  geom_label_repel(
    data = CD8.markers[which(CD8.markers$gene %in% genes.CD8), ], aes(x = avg_log2FC, y = FDR, label = gene),
    color = "black", label.size = 0, size = 2.5, max.overlaps = 20
  ) +
  scale_x_continuous("log2FC", limits = c(-2, 2)) +
  scale_y_continuous("-log10(FDR)", limits = c(0, 200)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/volcano/20241023_eAML_combined_CD8_BM_eAML.svg", width = 3, height = 3)


eAML.data.Tcells.CD4 <- subset(eAML.data.Tcells, cells = colnames(eAML.data.Tcells)[which(eAML.data.Tcells$manual.cluster == "CD4+ T cell")])
DefaultAssay(eAML.data.Tcells.CD4) <- "RNA"
eAML.data.Tcells.CD4 <- NormalizeData(eAML.data.Tcells.CD4)
eAML.data.Tcells.CD4 <- ScaleData(eAML.data.Tcells.CD4, features = rownames(eAML.data.Tcells.CD4))
eAML.data.Tcells.CD4 <- FindVariableFeatures(eAML.data.Tcells.CD4)
eAML.data.Tcells.CD4 <- RunPCA(eAML.data.Tcells.CD4)
eAML.data.Tcells.CD4 <- RunUMAP(eAML.data.Tcells.CD4, dims = 1:10)

CD4.markers <- FindMarkers(eAML.data.Tcells.CD4, group.by = "tissue", ident.1 = "BM", ident.2 = "eAML", logfc.threshold = 0.1)

genes.CD4 <- c(
  "ZFP36", "CCL4", "CXCR4", "ZFP36L1", "ZNF331", "CD44", "ZFP36L2", "CTLA4", "ICOS", "IFNG", "PDCD1",
  "ZAP70", "ITGA4", "ITGB2", "SELL", "BCL2", "TCF7", "ITGA6", "ITGB7"
)

CD4.markers$FDR <- -log10(CD4.markers$p_val_adj)
CD4.markers <- CD4.markers[which(abs(CD4.markers$avg_log2FC) < 2 & CD4.markers$FDR < 150), ]
CD4.markers$gene <- rownames(CD4.markers)

ggplot() +
  ggrastr::rasterize(geom_point(data = CD4.markers, aes(x = avg_log2FC, y = FDR), color = "grey", size = 0.5), dpi = 600) +
  ggrastr::rasterize(geom_point(data = CD4.markers[which(CD4.markers$gene %in% genes.CD4), ], aes(x = avg_log2FC, y = FDR), color = "firebrick", size = 0.5), dpi = 600) +
  geom_label_repel(
    data = CD4.markers[which(CD4.markers$gene %in% genes.CD4), ], aes(x = avg_log2FC, y = FDR, label = gene),
    color = "black", label.size = 0, size = 2.5, max.overlaps = 20
  ) +
  scale_x_continuous("log2FC", limits = c(-2, 2)) +
  scale_y_continuous("-log10(FDR)", limits = c(0, 150)) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/volcano/20241023_eAML_combined_CD4_BM_eAML.svg", width = 3, height = 3)

### genes of interest

genes <- c(
  "CD3E", "CD4", "CD8A", "SELL", "CCR7", "IL7R", "CD28", "FAS", "CD27", "ITGA4", "ITGAE", "ITGAL", "ITGAM", "ITGAX", "ITGB2",
  "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4", "CD226", "VTCN1", "CD244", "KLRG1", "TNFRSF14", "BTLA", "CD160",
  "CD38", "ENTPD1", "NT5E", "CD69", "IL2RA", "ICOS", "TNFRSF4", "TNFRSF9", "HLA-DRA", "CD40LG",
  "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "NKG7", "GNLY", "IFNG", "FASLG", "TNF", "IL17A", "IL2",
  "LEF1", "TCF7", "EOMES", "TBX21", "PRDM1", "TOX", "GATA3", "ID2", "ID3", "NR4A1", "ZNF683", "FOXP3"
)

expr.mat <- data.frame()

eAML.data.Tcells <- ScaleData(eAML.data.Tcells, features = genes)
expr.mat <- as.data.frame(t(
  GetAssayData(eAML.data.Tcells,
    layer = "scale.data"
  )[
    genes,
    colnames(eAML.data.Tcells)[which(eAML.data.Tcells$manual.cluster %in%
      c("CD4+ T cell", "CD8+ T cell", "Treg"))]
  ]
))
expr.mat$tissue <- eAML.data.Tcells$tissue[rownames(expr.mat)]
expr.mat$manual.cluster <- eAML.data.Tcells$manual.cluster[rownames(expr.mat)]

expr.mat <- expr.mat %>%
  group_by(tissue, manual.cluster) %>%
  summarize_all("mean")

ha <- columnAnnotation(
  tissue = expr.mat$tissue,
  col = list(
    "tissue" = c("BM" = "firebrick", "PB" = "red", "eAML" = "orange", "skin" = "blue"),
    "samples" = c(BM.samples.colors, eAML.samples.colors, PB.samples.colors, skin.samples.colors)
  ),
  border = T,
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 7), simple_anno_size = unit(5, "pt")
)



col_fun <- circlize::colorRamp2(breaks = seq(-2, 2, 4 / 8), colors = BuenColors::jdb_palette(name = "brewer_yes"))
svglite::svglite("./scRNA_seq/figures/combined/heatmap/20241023_expression.svg", width = 4, height = 5.5)
Heatmap(as.matrix(t(expr.mat[3:ncol(expr.mat)])), # column_split = factor(eAML.data.Tcells.CD4$orig.ident, levels = c(BM.samples, PB.samples, eAML.samples, skin.samples)),
  # column_split = factor(expr.mat$tissue, levels = c('BM', 'PB', 'eAML', 'skin')),
  column_split = factor(expr.mat$manual.cluster, levels = c("CD4+ T cell", "CD8+ T cell", "Treg")),
  cluster_rows = F, cluster_columns = F, column_title_gp = gpar(fontsize = 7), col = col_fun, row_names_gp = gpar(fontsize = 7),
  show_row_dend = F, show_column_dend = F, show_column_names = F, row_names_side = "left", row_title_gp = gpar(fontsize = 7),
  top_annotation = ha, use_raster = T, raster_quality = 10, border = T
)
dev.off()


#### exhaustion and memory score

DefaultAssay(eAML.data.Tcells) <- "RNA"
eAML.data.Tcells <- AddModuleScore(eAML.data.Tcells,
  features = list("Oliveira.tumor.score" = c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "TOX", "LAG3", "ENTPD1")),
  ctrl = 5, name = "Oliveira.tumor.score"
)
eAML.data.Tcells <- AddModuleScore(eAML.data.Tcells,
  features = list("Oliveira.memory.score" = c("TCF7", "IL7R", "SELL", "CCR7", "CD28")),
  ctrl = 5, name = "Oliveira.memory.score"
)

df <- data.frame(
  tumor.score = eAML.data.Tcells$Oliveira.tumor.score1,
  memory.score = eAML.data.Tcells$Oliveira.memory.score1,
  manual.cluster = eAML.data.Tcells$manual.cluster,
  tissue = eAML.data.Tcells$tissue
)

df$exhausted <- ifelse(df$tumor.score > 0.5 & df$memory.score < 0.5, "exhausted", "non.exhausted")

df.stats <- df %>%
  group_by(manual.cluster, tissue) %>%
  summarize(exhausted.freq = length(which(exhausted == "exhausted")) /
    length(exhausted))


p <- ggplot(df[which(df$tissue == "BM" & df$manual.cluster == "CD4+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_BM_CD4.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "PB" & df$manual.cluster == "CD4+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_PB_CD4.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "eAML" & df$manual.cluster == "CD4+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_eAML_CD4.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "skin" & df$manual.cluster == "CD4+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_skin_CD4.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))


p <- ggplot(df[which(df$tissue == "BM" & df$manual.cluster == "CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_BM_CD8.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "PB" & df$manual.cluster == "CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_PB_CD8.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "eAML" & df$manual.cluster == "CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_eAML_CD8.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "skin" & df$manual.cluster == "CD8+ T cell"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_skin_CD8.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))



p <- ggplot(df[which(df$tissue == "BM" & df$manual.cluster == "Treg"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_BM_Treg.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "PB" & df$manual.cluster == "Treg"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_PB_Treg.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "eAML" & df$manual.cluster == "Treg"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_eAML_Treg.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$tissue == "skin" & df$manual.cluster == "Treg"), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20241023_exhaustion_memory_skin_Treg.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))
