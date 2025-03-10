# analysis of progenitor compartment

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(Seurat)


source("./analysis/00_eAML.colors.R")
source("./analysis/00_genes_of_interest.R")

eAML.data <- readRDS("./data/scRNA_seq/objects/20241021_eAML_combined.scvi.rds")

# subset progenitors and reintegrate with scvi
eAML.data.progenitor <- subset(eAML.data, manual.cluster %in% c("myeloid"))

library(reticulate)
library(anndata)
library(sceasy)

# need to adjust to local environment
use_condaenv("/Users/liviuspenter/miniconda3/envs/scvi-env/")

sc <- import("scanpy")
scvi <- import("scvi")

eAML.data.progenitor[["RNA"]] <- as(object = eAML.data.progenitor[["RNA"]], Class = "Assay")
adata <- convertFormat(eAML.data.progenitor, from = "seurat", to = "anndata", main_layer = "counts", drop_single_values = FALSE)

scvi$model$SCVI$setup_anndata(adata, batch_key = "orig.ident")
model <- scvi$model$SCVI(adata)
model$train()
latent <- model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) <- colnames(eAML.data.progenitor)
eAML.data.progenitor[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(eAML.data.progenitor))
eAML.data.progenitor <- FindNeighbors(eAML.data.progenitor, reduction = "scvi")
eAML.data.progenitor <- FindClusters(eAML.data.progenitor, resolution = 0.1)
eAML.data.progenitor <- RunUMAP(eAML.data.progenitor, dims = 1:10, reduction = "scvi", n.components = 2)
saveRDS(eAML.data.progenitor, file = "./data/scRNA_seq/objects/20241024_eAML_combined_progenitor.rds")
eAML.data.progenitor <- readRDS(file = "./data/scRNA_seq/objects/20241024_eAML_combined_progenitor.rds")

# use reference annotation with label transfer to help with identification of clusters
bm <- readRDS(file = "./data/scRNA_seq/objects/bm.reference.rds")
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./data/scRNA_seq/objects/reftmp.idx")

anchors <- FindTransferAnchors(
  reference = bm, query = eAML.data.progenitor, k.filter = NA,
  reference.reduction = "spca", reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)

eAML.data.progenitor <- MapQuery(
  anchorset = anchors, query = eAML.data.progenitor,
  reference = bm,
  refdata = list(
    celltype = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# overview UMAPs
p <- DimPlot(subset(eAML.data.progenitor, tissue %in% c("BM", "eAML", "PB")), shuffle = T, group.by = "tissue", cols = tissue.colors) +
  NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241105_progenitor_BM_eAML.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(subset(eAML.data.progenitor, tissue %in% c("BM", "eAML", "PB")), shuffle = T, group.by = "chimerism") +
  NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241105_progenitor_BM_eAML.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(eAML.data.progenitor,
  cols.highlight = donor.recipient.colors, sizes.highlight = 2,
  cells.highlight =
    list(
      "recipient" = colnames(eAML.data.progenitor)[which(eAML.data.progenitor$tissue == "PB" & eAML.data.progenitor$patient == "1003" &
        eAML.data.progenitor$chimerism == "recipient")],
      "donor" = colnames(eAML.data.progenitor)[which(eAML.data.progenitor$tissue == "PB" & eAML.data.progenitor$patient == "1003" &
        eAML.data.progenitor$chimerism == "donor")]
    )
) +
  NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241105_eAML2_PB.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(eAML.data.progenitor,
  cols.highlight = donor.recipient.colors["donor"], sizes.highlight = 2,
  cells.highlight =
    list(
      "recipient" = colnames(eAML.data.progenitor)[which(eAML.data.progenitor$tissue == "PB" & eAML.data.progenitor$patient == "EKP" &
        eAML.data.progenitor$chimerism == "recipient")],
      "donor" = colnames(eAML.data.progenitor)[which(eAML.data.progenitor$tissue == "PB" & eAML.data.progenitor$patient == "EKP" &
        eAML.data.progenitor$chimerism == "donor")]
    )
) +
  NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241105_eAML4_PB.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(eAML.data.progenitor,
  cols.highlight = donor.recipient.colors, sizes.highlight = 2,
  cells.highlight =
    list(
      "recipient" = colnames(eAML.data.progenitor)[which(eAML.data.progenitor$tissue == "PB" & eAML.data.progenitor$patient == "KB" &
        eAML.data.progenitor$chimerism == "recipient")],
      "donor" = colnames(eAML.data.progenitor)[which(eAML.data.progenitor$tissue == "PB" & eAML.data.progenitor$patient == "KB" &
        eAML.data.progenitor$chimerism == "donor")]
    )
) +
  NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/combined/UMAP/20241105_eAML7_PB.png", width = 4, height = 4, dpi = 600, plot = p)


### DGEA between bone marrow and eAML
# module scores for HLA expression
eAML.data.progenitor <- AddModuleScore(eAML.data.progenitor,
  features = list("HLA.II" = HLA.II),
  ctrl = 5, name = "HLA.II"
)
eAML.data.progenitor <- AddModuleScore(eAML.data.progenitor,
  features = list("HLA.I" = HLA.I),
  ctrl = 5, name = "HLA.I"
)

# compare only cells in shared cluster
boo <- subset(eAML.data.progenitor, manual.cluster %in% c("progenitor1", "progenitor2") & tissue %in% c("BM", "eAML"))

markers <- FindMarkers(boo, group.by = "tissue", ident.1 = "BM", ident.2 = "eAML")
markers$gene <- rownames(markers)
markers$FDR <- -log10(markers$p_val_adj)
markers$FDR[which(markers$p_val_adj == 0)] <- -log10(min(markers$p_val_adj[which(markers$p_val_adj != 0)]))

p <- ggplot() +
  ggrastr::rasterize(geom_point(data = markers, aes(x = avg_log2FC, y = FDR), size = 0.5, color = "lightblue", pch = 16, alpha = 0.2), dpi = 600) +
  geom_point(data = markers[which(markers$gene %in% c(HLA.I)), ], aes(x = avg_log2FC, y = FDR), size = 0.5, color = "blue") +
  geom_point(data = markers[which(markers$gene %in% c(HLA.II)), ], aes(x = avg_log2FC, y = FDR), size = 0.5, color = "firebrick") +
  scale_x_continuous("log2(FC)", limits = c(-5, 5)) +
  scale_y_continuous("-log10(FDR)") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/volcano/20241103_eAML_BM.svg", width = 2.5, height = 2, plot = p)

### GO enrichment
library(org.Hs.eg.db)
library(goseq)

genes <- rownames(markers)[which(markers$avg_log2FC > 0.5)]
genes <- ifelse(rownames(eAML.data.progenitor) %in% genes, 1, 0)
names(genes) <- rownames(eAML.data.progenitor)
pwf <- nullp(genes, "hg38", "geneSymbol")
goResults <- goseq(pwf, "hg38", "geneSymbol")
goResults$count <- 100 * goResults$numDEInCat / goResults$numInCat
goResults$p.adj <- p.adjust(goResults$over_represented_pvalue)
ggplot(
  goResults[which(goResults$count > 30 & goResults$p.adj < 10^-2), ],
  aes(x = reorder(term, count), y = count)
) +
  geom_point(aes(size = numDEInCat, fill = log10(p.adj)), color = "black", pch = 21) +
  coord_flip() +
  scale_fill_gradientn(
    colours = rev(BuenColors::jdb_palette(name = "solar_rojos")),
    limits = c(min(log10(goResults[which(goResults$count > 30 & goResults$p.adj < 10^-2), "p.adj"])), 0)
  ) +
  scale_y_continuous("% hits") +
  scale_size_continuous(range = c(1, 3.5)) +
  theme_bw() +
  theme( # legend.position = 'None',
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, color = "black"),
    axis.text.y = element_text("Arial", size = 8, color = "black")
  )
ggsave("./figures/combined/GO_enrichment/20241103_progenitor.svg", width = 7, height = 3)

HLA.expr <- data.frame(
  cb = colnames(boo),
  HLA.I = boo$HLA.I1,
  HLA.II = boo$HLA.II1,
  tissue = boo$tissue,
  sample = boo$orig.ident,
  patient = boo$patient,
  chimerism = boo$chimerism
)

HLA.expr <- HLA.expr[-which(HLA.expr == "donor"), ]

p <- ggplot(HLA.expr[sample(nrow(HLA.expr), nrow(HLA.expr), replace = F), ], aes(x = HLA.I, y = HLA.II)) +
  scale_x_continuous("HLA class I", limits = c(-1, 2)) +
  scale_y_continuous("HLA class II", limits = c(-1, 2)) +
  ggrastr::rasterize(geom_point(aes(color = tissue), alpha = 0.05, stroke = 0, pch = 16), dpi = 600) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_density_2d(data = HLA.expr[which(HLA.expr$tissue == "BM"), ], size = 0.75, color = "firebrick") +
  geom_density_2d(data = HLA.expr[which(HLA.expr$tissue == "eAML"), ], size = 0.75, color = "orange") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_HLA_I_II.svg", width = 2.5, height = 2.5, plot = p)

HLA.expr <- HLA.expr %>%
  group_by(tissue, patient) %>%
  summarize(
    HLA.I = mean(HLA.I),
    HLA.II = mean(HLA.II)
  ) %>%
  as.data.frame()

HLA.expr <- HLA.expr[order(HLA.expr$HLA.II), ]
HLA.expr$patient <- factor(HLA.expr$patient, levels = HLA.expr$patient)

p <- ggplot() +
  geom_col(data = HLA.expr, aes(x = patient, y = HLA.II, fill = tissue), color = "black") +
  scale_fill_manual(values = tissue.colors) +
  scale_x_discrete("AML case") +
  scale_y_continuous("HLA class II") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20241113_HLA_II_expression_by_case.svg", width = 2.5, height = 1.5)


q <- ggplot() +
  geom_col(data = HLA.expr, aes(x = patient, y = HLA.I, fill = tissue), color = "black") +
  scale_fill_manual(values = tissue.colors) +
  scale_x_discrete("AML case") +
  scale_y_continuous("HLA class I") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20241113_HLA_I_expression_by_case.svg", width = 2.5, height = 1.5)

plots <- cowplot::plot_grid(plotlist = list(q, p), ncol = 1, align = "hv")
ggsave("./scRNA_seq/figures/combined/plots/20241113_HLA_I_II_expression_by_case.svg", width = 2.5, height = 2, plot = plots)

# homing receptors
boo <- subset(eAML.data.progenitor, manual.cluster %in% c("progenitor1", "progenitor2") & tissue %in% c("BM", "eAML"))

boo <- ScaleData(boo, features = c(checkpoint.molecules, integrins, adhesion, selectins))
df <- GetAssayData(boo, layer = "scale.data")[c(checkpoint.molecules, integrins, adhesion, selectins), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame()
df$HLA.I <- boo$HLA.I1
df$HLA.II <- boo$HLA.II1
df$patient <- as.character(boo$patient)
df$tissue <- as.character(boo$tissue)
df <- df %>%
  group_by(patient, tissue) %>%
  summarize_all(.funs = mean)

ha <- columnAnnotation(tissue = df$tissue)

col_fun <- circlize::colorRamp2(breaks = seq(-4, 4, 8 / 8), colors = BuenColors::jdb_palette(name = "brewer_yes"))
Heatmap(scale(t(df[, which(!colnames(df) %in% c("patient", "tissue"))])),
  top_annotation = ha, column_split = df$tissue,
  row_order = c("HLA.I", "HLA.II", c(checkpoint.molecules, integrins, adhesion, selectins)), col = col_fun
)


ggplot(df, aes(x = tissue, y = ICAM1, color = tissue)) +
  geom_violin(width = 0.7, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_ICAM1.svg", width = 1.2, height = 2)

ggplot(df, aes(x = tissue, y = PECAM1, color = tissue)) +
  geom_violin(width = 0.7, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_PECAM1.svg", width = 1.2, height = 2)

ggplot(df, aes(x = tissue, y = ITGAL, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_ITGAL.svg", width = 1.3, height = 2)

ggplot(df, aes(x = tissue, y = ITGA4, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_ITGA4.svg", width = 1.3, height = 2)

ggplot(df, aes(x = tissue, y = ITGA5, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_ITGA5.svg", width = 1.3, height = 2)


ggplot(df, aes(x = tissue, y = ITGA6, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_ITGA6.svg", width = 1.3, height = 2)

ggplot(df, aes(x = tissue, y = ITGAV, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_ITGAV.svg", width = 1.3, height = 2)


ggplot(df, aes(x = tissue, y = ITGB4, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_ITGB4.svg", width = 1.3, height = 2)




ggplot(df, aes(x = tissue, y = LGALS9, color = tissue)) +
  geom_violin(width = 1) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.1, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_LGALS9.svg", width = 1.3, height = 2)

ggplot(df, aes(x = tissue, y = HMGB1, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_HMGB1.svg", width = 1.3, height = 2)

ggplot(df, aes(x = tissue, y = CEACAM1, color = tissue)) +
  geom_violin(width = 1, aes(fill = tissue), color = "black") +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "black") +
  geom_jitter(color = "black", width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(c("eAML", "BM")), color = "black", test = "t.test") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/combined/plots/20241103_CEACAM1.svg", width = 1.3, height = 2)


boo <- subset(eAML.data.progenitor, patient %in% c("1003", "EKP", "KB"))
boo <- ScaleData(boo, features = c(checkpoint.molecules, integrins, adhesion, selectins))
df <- GetAssayData(boo, layer = "scale.data")[c(checkpoint.molecules, integrins, adhesion, selectins), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame()
df$HLA.I <- boo$HLA.I1
df$HLA.II <- boo$HLA.II1
df$patient <- as.character(boo$patient)
df$tissue <- as.character(boo$tissue)
df$chimerism <- as.character(boo$chimerism)
df$patient.tissue.chimerism <- paste0(df$patient, ".", df$tissue, ".", df$chimerism)

mat <- df[which(df$chimerism == "recipient"), c("ITGA4", "ITGA5", "ITGA6", "ITGAL", "ITGAV", "ITGB4", "PECAM1", "ICAM1", "tissue")]

mat <- df[which(df$chimerism == "recipient"), c("ITGA4", "ITGA6", "ITGAL", "ITGAV", "PECAM1", "ICAM1", "tissue")]

ha <- columnAnnotation(
  tissue = mat$tissue, col = list("tissue" = tissue.colors), border = T,
  annotation_name_side = "left", simple_anno_size = unit(5, "pt")
)
col_fun <- circlize::colorRamp2(breaks = seq(-2, 2, 4 / 8), colors = BuenColors::jdb_palette(name = "brewer_yes"))
svglite::svglite("./scRNA_seq/figures/combined/heatmap/20241113_eAML_PB.svg", width = 3.5, height = 1.5)
Heatmap(t(mat[, 1:6]),
  column_split = factor(mat$tissue, levels = c("PB", "eAML")), show_row_dend = F, show_column_dend = F,
  cluster_columns = F, show_column_names = F, cluster_rows = T, border = T, col = col_fun, top_annotation = ha,
  row_names_side = "left", raster_quality = 10, use_raster = T
)
dev.off()

boo$tissue <- factor(boo$tissue, levels = c("PB", "eAML"))



for (g in c("ITGA4", "ITGA5", "ITGA6", "ITGAL", "ITGAV", "ITGB4", "PECAM1", "ICAM1")) {
  VlnPlot(subset(boo, chimerism == "recipient" & patient %in% c("1003", "KB")), group.by = "patient", split.by = "tissue", features = g) +
    scale_fill_manual(values = tissue.colors) + NoLegend() +
    scale_x_discrete(labels = c("eAML2", "eAML7")) +
    scale_y_continuous(paste0(g)) +
    theme(
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_blank()
    )
  ggsave(paste0("./scRNA_seq/figures/combined/plots/20241216_", g, "_PB_BM.svg"), width = 1.3, height = 2)
}
