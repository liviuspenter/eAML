# QC metrics for comparison of different tissue-sites

source("./scRNA_seq/R/00_eAML.colors.R")

library(dplyr)
library(ggplot2)
library(ggsignif)
library(Seurat)

eAML.combined <- readRDS("./data/scRNA_seq/objects/20241021_eAML_combined.scvi.rds")

unique(eAML.combined$orig.ident)

df <- data.frame(
  cb = colnames(eAML.combined),
  percent.mt = eAML.combined$percent.mt,
  nCount_RNA = eAML.combined$nCount_RNA,
  nFeature_RNA = eAML.combined$nFeature_RNA,
  orig.ident = eAML.combined$orig.ident,
  tissue = eAML.combined$tissue
)

stats <- df %>%
  group_by(orig.ident, tissue) %>%
  summarize(
    percent.mt = mean(percent.mt),
    nCount_RNA = mean(nCount_RNA),
    nFeature_RNA = mean(nFeature_RNA),
    cells = length(orig.ident)
  )

stats$tissue <- factor(stats$tissue, levels = c("BM", "PB", "eAML", "skin"))

ggplot(stats, aes(x = tissue, y = percent.mt, fill = tissue)) +
  geom_boxplot(outlier.size = 0, color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_manual(values = tissue.colors) +
  scale_y_continuous("%mitochondrial transcripts") +
  geom_signif(comparisons = list(c("BM", "PB"), c("BM", "eAML"), c("BM", "skin"), c("eAML", "skin")), step_increase = 0.15, textsize = 3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/qc/20250414_eAML_combined_percent_mt.svg", width = 1.3, height = 2.5)

ggplot(stats, aes(x = tissue, y = nCount_RNA, fill = tissue)) +
  geom_boxplot(outlier.size = 0, color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_manual(values = tissue.colors) +
  scale_y_continuous("UMI per cell") +
  geom_signif(comparisons = list(c("BM", "PB"), c("BM", "eAML"), c("BM", "skin"), c("eAML", "skin")), step_increase = 0.15, textsize = 3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/qc/20250414_eAML_combined_nCount_RNA.svg", width = 1.5, height = 2.5)

ggplot(stats, aes(x = tissue, y = nFeature_RNA, fill = tissue)) +
  geom_boxplot(outlier.size = 0, color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_manual(values = tissue.colors) +
  scale_y_continuous("genes per cell") +
  geom_signif(comparisons = list(c("BM", "PB"), c("BM", "eAML"), c("BM", "skin"), c("eAML", "skin")), step_increase = 0.15, textsize = 3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/qc/20250414_eAML_combined_nFeature_RNA.svg", width = 1.5, height = 2.5)

ggplot(stats, aes(x = tissue, y = cells, fill = tissue)) +
  geom_boxplot(outlier.size = 0, color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_fill_manual(values = tissue.colors) +
  scale_y_continuous("cells per sample") +
  geom_signif(comparisons = list(c("BM", "PB"), c("BM", "eAML"), c("BM", "skin"), c("eAML", "skin")), step_increase = 0.15, textsize = 3) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./scRNA_seq/figures/qc/20250414_eAML_combined_cells.svg", width = 1.5, height = 2.5)
