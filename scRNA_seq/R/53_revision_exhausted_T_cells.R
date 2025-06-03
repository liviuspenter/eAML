# absolute number of exhausted T cells across tissues

library(dplyr)
library(ggplot2)
library(ggsignif)
library(Seurat)
library(tidyr)

source("./scRNA_seq/R/00_eAML.colors.R")

eAML.data <- readRDS("./data/scRNA_seq/objects/20241021_eAML_combined.scvi.rds")
eAML.data.Tcells <- readRDS("./data/scRNA_seq/objects/20241017_eAML_combined_Tcells.rds")

# calculate T cell exhaustion
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
  sample = eAML.data.Tcells$orig.ident,
  tissue = eAML.data.Tcells$tissue
)

df$exhausted <- ifelse(df$tumor.score > 0.5 & df$memory.score < 0.5, "exhausted", "non.exhausted")

# df.stats = df %>% group_by(manual.cluster, sample) %>%
#  summarize(exhausted.freq = length(which(exhausted == 'exhausted')) /
#              length(exhausted))
# df.stats = df.stats %>% pivot_wider(names_from = manual.cluster, values_from = exhausted.freq)

df.all <- data.frame(
  bc = colnames(eAML.data),
  sample = eAML.data$orig.ident,
  manual.cluster = eAML.data$manual.cluster,
  tissue = eAML.data$tissue
)

stats.all <- df.all %>%
  # filter(!sample %in% c('SKN8104902', 'SKN8105200')) %>%
  group_by(sample, manual.cluster, tissue) %>%
  summarize(cells = length(sample)) %>%
  pivot_wider(names_from = manual.cluster, values_from = cells)
stats.all$CD4 <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD4+ T cell"))
})
stats.all$CD8 <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD8+ T cell"))
})
stats.all$Treg <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "Treg"))
})
stats.all$exhausted.CD4 <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD4+ T cell" & df$exhausted == "exhausted"))
})
stats.all$exhausted.CD8 <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD8+ T cell" & df$exhausted == "exhausted"))
})
stats.all$exhausted.Treg <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "Treg" & df$exhausted == "exhausted"))
})


stats.all$exhausted.CD4.freq <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD4+ T cell" & df$exhausted == "exhausted")) /
    length(which(df$sample == x))
})
stats.all$exhausted.CD8.freq <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD8+ T cell" & df$exhausted == "exhausted")) /
    length(which(df$sample == x))
})
stats.all$exhausted.Treg.freq <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "Treg" & df$exhausted == "exhausted")) /
    length(which(df$sample == x))
})

stats.all$CD4.freq <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD4+ T cell")) /
    length(which(df$sample == x))
})
stats.all$CD8.freq <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "CD8+ T cell")) /
    length(which(df$sample == x))
})
stats.all$Treg.freq <- sapply(stats.all$sample, FUN = function(x) {
  length(which(df$sample == x & df$manual.cluster == "Treg")) /
    length(which(df$sample == x))
})


stats.all$tissue <- factor(stats.all$tissue, levels = c("BM", "PB", "eAML", "skin"))
stats.all$myeloid.exhausted.CD4 <- stats.all$myeloid / stats.all$exhausted.CD4
stats.all$myeloid.exhausted.CD8 <- stats.all$myeloid / stats.all$exhausted.CD8
stats.all$myeloid.exhausted.Treg <- stats.all$myeloid / stats.all$exhausted.Treg
stats.all$myeloid.CD4 <- stats.all$myeloid / stats.all$CD4
stats.all$myeloid.CD8 <- stats.all$myeloid / stats.all$CD8
stats.all$myeloid.Treg <- stats.all$myeloid / stats.all$Treg

ggplot(stats.all[which(stats.all$myeloid.CD4 != Inf), ], aes(x = tissue, y = myeloid.CD4)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("myeloid per CD4+ T cell") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_myeloid_CD4_ratio.svg", width = 1.5, height = 2.5)

ggplot(stats.all[which(stats.all$myeloid.exhausted.CD4 != Inf), ], aes(x = tissue, y = myeloid.exhausted.CD4)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("myeloid per exh. CD4+ T cell") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_myeloid_exhausted_CD4_ratio.svg", width = 1.5, height = 2.5)

ggplot(stats.all[which(!is.na(stats.all$exhausted.CD4.freq)), ], aes(x = tissue, y = 100 * exhausted.CD4.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("% exh. CD4+ T cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_exhausted_CD4_freq.svg", width = 1.5, height = 2.5)

ggplot(stats.all[which(!is.na(stats.all$CD4.freq)), ], aes(x = tissue, y = 100 * CD4.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("% CD4+ T cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_CD4_freq.svg", width = 1.5, height = 2.5)


ggplot(stats.all[which(stats.all$myeloid.exhausted.CD8 != Inf), ], aes(x = tissue, y = myeloid.exhausted.CD8)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("myeloid per exh. CD8+ T cell") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_myeloid_exhausted_CD8_ratio.svg", width = 1.5, height = 2.5)

ggplot(stats.all[which(stats.all$myeloid.CD8 != Inf), ], aes(x = tissue, y = myeloid.CD8)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("myeloid per CD8+ T cell") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_myeloid_CD8_ratio.svg", width = 1.5, height = 2.5)

ggplot(stats.all[which(!is.na(stats.all$exhausted.CD8.freq)), ], aes(x = tissue, y = 100 * exhausted.CD8.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("% exh. CD8+ T cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_exhausted_CD8_freq.svg", width = 1.5, height = 2.5)

ggplot(stats.all[which(!is.na(stats.all$CD8.freq)), ], aes(x = tissue, y = 100 * CD8.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("% CD8+ T cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_CD8_freq.svg", width = 1.5, height = 2.5)


ggplot(stats.all[which(stats.all$myeloid.exhausted.Treg != Inf), ], aes(x = tissue, y = myeloid.exhausted.Treg)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("myeloid per exh. Treg") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_myeloid_exhausted_Treg_ratio.svg", width = 1.5, height = 2.5)


ggplot(stats.all[which(stats.all$myeloid.Treg != Inf), ], aes(x = tissue, y = myeloid.Treg)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("myeloid per Treg") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_myeloid_Treg_ratio.svg", width = 1.5, height = 2.5)



ggplot(stats.all[which(!is.na(stats.all$exhausted.Treg.freq)), ], aes(x = tissue, y = 100 * exhausted.Treg.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("% exh. Tregs") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_exhausted_Treg_freq.svg", width = 1.5, height = 2.5)

ggplot(stats.all[which(!is.na(stats.all$Treg.freq)), ], aes(x = tissue, y = 100 * Treg.freq)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1) +
  geom_signif(comparisons = list(c("BM", "eAML"), c("eAML", "skin"), c("BM", "skin")), step_increase = 0.15, textsize = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_y_continuous("% Tregs") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/combined/plots/20250423_Treg_freq.svg", width = 1.5, height = 2.5)
