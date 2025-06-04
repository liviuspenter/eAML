# compare T cell phenotypes between bone marrow and extramedullary AML in eAML1

library(circlize)
library(colorRamps)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(Seurat)

source("./scRNA_seq/R/00_eAML.colors.R")

genes <- c(
  "CD3E", "CD4", "CD8A", "SELL", "CCR7", "IL7R", "CD28", "FAS", "CD27", "ITGA4", "ITGA5", "ITAG6", "ITGAE", "ITGAL", "ITGAM", "ITGAV5", "ITGAX", "ITGB2", "ITGB4",
  "PECAM1", "ICAM1",
  "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4", "CD226", "VTCN1", "CD244", "KLRG1", "TNFRSF14", "BTLA", "CD160",
  "CD38", "ENTPD1", "NT5E", "CD69", "IL2RA", "ICOS", "TNFRSF4", "TNFRSF9", "HLA-DRA", "CD40LG",
  "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "NKG7", "GNLY", "IFNG", "FASLG", "TNF", "IL17A", "IL2",
  "LEF1", "TCF7", "EOMES", "TBX21", "PRDM1", "TOX", "GATA3", "ID2", "ID3", "NR4A1", "ZNF683", "FOXP3"
)


AML44 <- readRDS("./data/scRNA_seq/objects/20250429_eAML1.combined.rds")
genes <- intersect(genes, rownames(AML44))
TCR.clones <- data.table::fread("./data/scRNA_seq/objects/20250429_eAML1.combined.TCR.clones.csv") %>% as.data.frame()


### circos plot between eAML1.1 and Pool97_23

choose_colour_clone <- function(x) {
  output <- c()
  for (i in 1:length(x)) {
    if (x[i] > 0) {
      output <- c(output, c(rep(primary.colors(24), ceiling(max(x) / 24)))[x[i]])
    }
    if (x[i] == 0) {
      output <- c(output, "grey")
    }
  }
  return(output)
}

# Define sectors (one for each sample)
sectors <- c("eAML1", "Pool97_23")

TCR.clones.reduced <- TCR.clones[which(TCR.clones$eAML1.1 != 0 | TCR.clones$Pool97_23 != 0), ]

# Normalize the cell counts into ranges for plotting
# Accumulate positions for the clones in each sector

TCR.clones.reduced <- TCR.clones.reduced[order(TCR.clones.reduced$eAML1.1, decreasing = T), ]
eAML1_positions <- data.frame(
  sector = "eAML1",
  start = cumsum(c(0, head(TCR.clones.reduced$eAML1.1, -1))),
  end = cumsum(TCR.clones.reduced$eAML1.1),
  clone = TCR.clones.reduced$CloneID
)

TCR.clones.reduced <- TCR.clones.reduced[order(TCR.clones.reduced$Pool97_23, decreasing = T), ]
Pool97_positions <- data.frame(
  sector = "Pool97_23",
  start = cumsum(c(0, head(TCR.clones.reduced$Pool97_23, -1))),
  end = cumsum(TCR.clones.reduced$Pool97_23),
  clone = TCR.clones.reduced$CloneID
)

# Combine all positions
all_pos <- rbind(eAML1_positions, Pool97_positions)

circos.data <- data.frame(
  well.number = c(
    seq(1, max(eAML1_positions$end)),
    seq(1, max(Pool97_positions$end))
  ),
  sample = c(
    rep("eAML1", max(eAML1_positions$end)),
    rep("Pool97_23", max(Pool97_positions$end))
  )
)
circos.data$tissue.color <- ifelse(circos.data$sample == "eAML1", tissue.colors["eAML"], tissue.colors["BM"])
circos.data$color <- NA
circos.data$clone <- NA

for (clone in TCR.clones.reduced[order(TCR.clones.reduced$eAML1.1, decreasing = T), "CloneID"]) {
  circos.data[
    which(circos.data$sample == "eAML1" & circos.data$well.number %in% seq(
      eAML1_positions$start[which(eAML1_positions$clone == clone)],
      eAML1_positions$end[which(eAML1_positions$clone == clone)]
    )),
    "color"
  ] <- ifelse(TCR.clones$eAML1.1[which(TCR.clones$CloneID == clone)] < 2, "grey", choose_colour_clone(clone))
  circos.data[
    which(circos.data$sample == "eAML1" & circos.data$well.number %in% seq(
      eAML1_positions$start[which(eAML1_positions$clone == clone)],
      eAML1_positions$end[which(eAML1_positions$clone == clone)]
    )),
    "clone"
  ] <- clone
}

for (clone in TCR.clones.reduced[order(TCR.clones.reduced$Pool97_23, decreasing = T), "CloneID"]) {
  circos.data[
    which(circos.data$sample == "Pool97_23" & circos.data$well.number %in% seq(
      Pool97_positions$start[which(Pool97_positions$clone == clone)],
      Pool97_positions$end[which(Pool97_positions$clone == clone)]
    )),
    "color"
  ] <- ifelse(TCR.clones$Pool97_23[which(TCR.clones$CloneID == clone)] < 2, "grey", choose_colour_clone(clone))
  circos.data[
    which(circos.data$sample == "Pool97_23" & circos.data$well.number %in% seq(
      Pool97_positions$start[which(Pool97_positions$clone == clone)],
      Pool97_positions$end[which(Pool97_positions$clone == clone)]
    )),
    "clone"
  ] <- clone
}
sample.levels <- factor(circos.data$sample, levels = c("eAML1", "Pool97_23"))

svglite::svglite("./scRNA_seq/figures/combined/circos/20250430_eAML1_circosplot.svg", width = 4, height = 4)
circos.par("track.height" = 0.05, start.degree = 90, gap.degree = 5)
circos.initialize(factors = sample.levels, x = circos.data$well.number)

circos.track(factors = sample.levels, y = circos.data$well.number, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), CELL_META$sector.index, cex = 0.5)
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  breaks <- seq(xlim[1], xlim[2], by = 1)
  n_breaks <- length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
    breaks[-1], rep(ylim[2], n_breaks - 1),
    col = circos.data$tissue.color[which(circos.data$sample == CELL_META$sector.index)], border = NA
  )
})

circos.par("track.height" = 0.1)

circos.track(factors = sample.levels, y = circos.data$well.number, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), CELL_META$sector.index, cex = 1.5)
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim
  breaks <- seq(xlim[1], xlim[2], by = 1)
  n_breaks <- length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
    breaks[-1], rep(ylim[2], n_breaks - 1),
    col = circos.data$color[which(circos.data$sample == CELL_META$sector.index)], border = NA
  )
})


for (clone in TCR.clones.reduced[order(TCR.clones.reduced$Pool97_23, decreasing = F), "CloneID"]) {
  if (TCR.clones.reduced$Pool97_23[which(TCR.clones.reduced$CloneID == clone)] != 0 &
    TCR.clones.reduced$eAML1.1[which(TCR.clones.reduced$CloneID == clone)] != 0) {
    message(clone)
    clone.color <- ifelse(TCR.clones$eAML1.1[which(TCR.clones$CloneID == clone)] > 1 | TCR.clones$Pool97_23[which(TCR.clones$CloneID == clone)] > 1, choose_colour_clone(clone), "grey")
    message(clone.color)

    circos.link("eAML1", c(
      eAML1_positions$start[which(eAML1_positions$clone == clone)],
      eAML1_positions$end[which(eAML1_positions$clone == clone)]
    ),
    "Pool97_23", c(
      Pool97_positions$start[which(Pool97_positions$clone == clone)],
      Pool97_positions$end[which(Pool97_positions$clone == clone)]
    ),
    col = clone.color, lwd = 0.01
    )
  }
}
dev.off()

### phenotype of shared T cell clones

shared.clones <- TCR.clones$CloneID[which(TCR.clones$eAML1.1 > 9 & TCR.clones$Pool97_23 != 0)]
shared.clones.colors <- BuenColors::jdb_palette(name = "corona", n = 9)

eAML.clones <- TCR.clones$CloneID[which(TCR.clones$eAML1.1 > 9 & TCR.clones$Pool97_23 == 0)]
BM.clones <- TCR.clones$CloneID[which(TCR.clones$eAML1.1 == 0 & TCR.clones$Pool97_23 > 10)]

names(shared.clones.colors) <- shared.clones
AML44 <- ScaleData(AML44, features = genes)
# gather data
expr.data <- data.frame()
for (clone in shared.clones) {
  boo <- rowMeans(GetAssayData(AML44, slot = "scale.data")[genes, colnames(AML44)[which(AML44$CloneID == clone & AML44$orig.ident == "eAML1.1")]]) %>%
    t() %>%
    as.data.frame()
  boo$CloneID <- clone
  boo$tissue <- "eAML"
  expr.data <- rbind(expr.data, boo)

  boo <- rowMeans(GetAssayData(AML44, slot = "scale.data")[genes, colnames(AML44)[which(AML44$CloneID == clone & AML44$orig.ident %in% c("Pool97_23", "Pool97_24", "Pool97_26"))]]) %>%
    t() %>%
    as.data.frame()
  boo$CloneID <- clone
  boo$tissue <- "BM"
  expr.data <- rbind(expr.data, boo)
}

for (gene in genes) {
  message(gene)
  p <- ggplot(expr.data, aes(x = tissue, y = !!sym(gene))) +
    geom_boxplot(outliers = F, aes(fill = tissue)) +
    geom_jitter(width = 0.2, size = 0.5) +
    geom_signif(comparisons = list(c("BM", "eAML")), textsize = 3) +
    scale_fill_manual(values = tissue.colors) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text("Arial", size = 10, color = "black"),
      axis.title = element_text("Arial", size = 10, color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  ggsave(paste0("./scRNA_seq/figures/eAML1/overlap_genes/20250513_", gene, ".svg"), width = 1.1, height = 1.7, plot = p)
}

for (clone in eAML.clones) {
  boo <- rowMeans(GetAssayData(AML44, slot = "scale.data")[genes, colnames(AML44)[which(AML44$CloneID == clone & AML44$orig.ident == "eAML1.1")]]) %>%
    t() %>%
    as.data.frame()
  boo$CloneID <- clone
  boo$tissue <- "eAML.only"
  expr.data <- rbind(expr.data, boo)
}

for (clone in BM.clones) {
  if (length(which(AML44$CloneID == clone)) != 0) {
    boo <- rowMeans(GetAssayData(AML44, slot = "scale.data")[genes, colnames(AML44)[which(AML44$CloneID == clone & AML44$orig.ident %in% c("Pool97_23", "Pool97_24", "Pool97_26"))]]) %>%
      t() %>%
      as.data.frame()
    boo$CloneID <- clone
    boo$tissue <- "BM.only"
    expr.data <- rbind(expr.data, boo)
  }
}

clone.colors <- as.character(BuenColors::jdb_palette(name = "corona", n = 21))
names(clone.colors) <- as.character(unique(expr.data$CloneID))

col_fun <- circlize::colorRamp2(breaks = seq(-2, 2, 4 / 8), colors = BuenColors::jdb_palette(name = "brewer_yes"))
ha <- columnAnnotation(
  tissue = expr.data$tissue,
  clone = expr.data$CloneID,
  col = list(
    "tissue" = c(tissue.colors, "BM.only" = as.character(tissue.colors["BM"]), "eAML.only" = as.character(tissue.colors["eAML"])),
    "clone" = clone.colors
  ),
  border = T, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10),
  annotation_label = c("tissue", "TCR clone"), simple_anno_size = unit(7, "pt")
)
mat <- t(expr.data[, 1:length(genes)])
svglite::svglite("./scRNA_seq/figures/combined/heatmap/20250429_comparison_Tcell_clones.svg", width = 5, height = 6.5)
Heatmap(mat,
  cluster_rows = F, cluster_columns = F, top_annotation = ha, column_split = factor(expr.data$tissue, levels = c("BM.only", "BM", "eAML", "eAML.only")),
  row_names_side = "left", col = col_fun, row_names_gp = gpar(fontsize = 8, fontface = "italic"), border = T
)
dev.off()


# phenotype of individual T cell clones


DefaultAssay(AML44) <- "RNA"
AML44 <- AddModuleScore(AML44,
  features = list("Oliveira.tumor.score" = c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "TOX", "LAG3", "ENTPD1")),
  ctrl = 5, name = "Oliveira.tumor.score"
)
AML44 <- AddModuleScore(AML44,
  features = list("Oliveira.memory.score" = c("TCF7", "IL7R", "SELL", "CCR7", "CD28")),
  ctrl = 5, name = "Oliveira.memory.score"
)

shared.clones <- TCR.clones.reduced$CloneID[which(TCR.clones.reduced$eAML1.1 != 0 & TCR.clones.reduced$Pool97_23 != 0)]

df <- data.frame(
  tumor.score = AML44$Oliveira.tumor.score1,
  memory.score = AML44$Oliveira.memory.score1,
  CloneID = AML44$CloneID,
  sample = AML44$orig.ident
)
df <- df[which(!is.na(df$CloneID)), ]
df$exhausted <- ifelse(df$tumor.score > 0.5 & df$memory.score < 0.5, "exhausted", "non.exhausted")
df$shared.clone <- ifelse(df$CloneID %in% shared.clones, "yes", "no")

df.stats <- df %>%
  group_by(sample, shared.clone) %>%
  summarize(exhausted.freq = length(which(exhausted == "exhausted")) /
    length(exhausted))

p <- ggplot(df[which(df$sample == "eAML1.1" & df$CloneID %in% shared.clones), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20250430_exhaustion_memory_eAML1.1_shared.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))

p <- ggplot(df[which(df$sample == "eAML1.1" & !df$CloneID %in% shared.clones), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20250430_exhaustion_memory_eAML1.1_unique.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))


p <- ggplot(df[which(df$sample == "Pool97_23" & df$CloneID %in% shared.clones), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20250430_exhaustion_memory_Pool97_23_shared.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))


p <- ggplot(df[which(df$sample == "Pool97_23" & !df$CloneID %in% shared.clones), ], aes(x = memory.score, y = tumor.score)) +
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
ggsave("./scRNA_seq/figures/combined/plots/20250430_exhaustion_memory_Pool97_23_unique.svg", width = 1.3, height = 1.3, plot = ggrastr::rasterize(p, dpi = 600))
