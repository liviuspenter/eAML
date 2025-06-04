# scRNA-seq analysis of eAML1 including TCR repertoires

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(Seurat)

source("./scRNA_seq/R/00_eAML.colors.R")

metadata <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 3))
rownames(metadata) <- metadata$sample
metadata <- metadata[which(metadata$patient == "44"), ]

eAML.genotype <- data.frame()
for (s in metadata$sample) {
  f <- as.data.frame(data.table::fread(paste0(metadata[s, "library"], "/souporcell/clusters.tsv")))
  f$barcode <- paste0(s, "_", f$barcode)
  eAML.genotype <- rbind(eAML.genotype, f)
}
rownames(eAML.genotype) <- eAML.genotype$barcode

metadata <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 1))
rownames(metadata) <- metadata$sample
metadata <- metadata[which(metadata$patient == "44"), ]

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

# exclude doublets
eAML.data <- subset(eAML.data, cells = eAML.genotype$barcode[which(eAML.genotype$status == "singlet")])

# read TCR
source("./scRNA_seq/R/00_read_clonotypes.R")
source("./scRNA_seq/R/00_read_tcr_bcr.R")

metadata <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 2))
rownames(metadata) <- metadata$sample
metadata <- metadata[which(metadata$patient == "44"), ]

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
eAML.data <- FindClusters(eAML.data, resolution = 0.3)
saveRDS("./data/scRNA_seq/objects/20230117_eAML1.rds", object = eAML.data)
eAML.data <- readRDS(file = "./data/scRNA_seq/objects/20230117_eAML1.rds")

p <- DimPlot(eAML.data, cols = eAML1.colors) +
  NoLegend() +
  NoAxes()
ggsave(filename = "./scRNA_seq/figures/eAML1/20241125_eAML1_UMAP.png", width = 4, height = 4, dpi = 600, plot = p)


eAML.data.markers <- FindAllMarkers(eAML.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(eAML.data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC))

eAML.data$manual.cluster <- "AML"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(1))] <- "T cell"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(5))] <- "fibroblast"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(6, 10, 12))] <- "keratinocyte"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(8))] <- "macrophage"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(9))] <- "endothelium"
eAML.data$manual.cluster[which(eAML.data$seurat_clusters %in% c(11))] <- "B cell"
eAML.data$manual.cluster <- factor(eAML.data$manual.cluster, levels = names(eAML1.colors))
DimPlot(eAML.data, group.by = "manual.cluster", cols = eAML1.colors) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./figures/eAML1/20230119_eAML1_celltypes.png", width = 4, height = 4, dpi = 600)

Idents(eAML.data) <- "manual.cluster"
df <- as.data.frame(prop.table(table(Idents(eAML.data), eAML.data$orig.ident), margin = 2))

ggplot(df, aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col(color = "black") +
  scale_x_discrete(labels = c("Pre-Nivo", "Post-Nivo", "Post-Ipi 1", "Post-Ipi 2", "inflammation", "Progression")) +
  scale_y_continuous("% cells") +
  scale_fill_manual("celltype", values = eAML1.colors) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_celltypes_kinetics.svg", width = 3, height = 2)

eAML.data$manual.cluster.chimerism <- paste0(eAML.data$manual.cluster, ".", eAML.data$chimerism)
eAML.data$manual.cluster.chimerism <- factor(eAML.data$manual.cluster.chimerism,
  levels = paste0(rep(names(eAML1.colors), each = 2), c(".donor", ".recipient"))
)
Idents(eAML.data) <- "manual.cluster.chimerism"
df <- as.data.frame(prop.table(table(Idents(eAML.data), eAML.data$orig.ident), margin = 2))

ggplot(df, aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col() +
  scale_x_discrete(labels = c("Pre-Nivo", "Post-Nivo", "Post-Ipi 1", "Post-Ipi 2", "inflammation", "Progression")) +
  scale_y_continuous("% cells") +
  scale_fill_manual("celltype", values = eAML1.colors.chimerism) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_celltypes_chimerism_kinetics.svg", width = 1.7, height = 2)

p <- DimPlot(eAML.data, group.by = "chimerism", cols = donor.recipient.colors) +
  NoLegend() +
  NoAxes() + theme(plot.title = element_blank())
ggsave(filename = "./scRNA_seq/figures/eAML1/20241125_eAML1_UMAP_chimerism.png", width = 4, height = 4, dpi = 600, plot = p)


# read TCR repertoire data and compare to scTCR-seq data
tcr.repertoire.data <- data.frame()
for (tcr_chain in c("alpha", "beta")) {
  for (s in c("44A", "44B", "44C", "44D", "44E", "44F", "44G")) {
    f <- read.csv(paste0("./data/bulk/TCRseq/tcr_", tcr_chain, "_repertoire_", s, ".csv"), row.names = 1)
    f$tcr_chain <- tcr_chain
    tcr.repertoire.data <- rbind(tcr.repertoire.data, f)
  }
}
tcr.repertoire.data$sample[which(tcr.repertoire.data$sample == "44")] <- "44E"
tcr.repertoire <- tcr.repertoire.data
tcr.repertoire.data$sample[which(grepl("44", tcr.repertoire.data$sample))] <- "eAML1"
tcr.repertoire.data <- as.data.frame(tcr.repertoire.data %>%
  group_by(sample, tcr_chain, CDR3) %>%
  summarize(UMI = sum(UMI)))

clonotypes$overlap <- ifelse(clonotypes$cdr3 %in% tcr.repertoire.data$CDR3, "yes", "no")
CDR3.overlap <- clonotypes %>%
  filter(cdr3 != "") %>%
  group_by(cdr3, chain) %>%
  summarize(n = n()) %>%
  group_by(chain) %>%
  mutate(rank = rank(-n, ties.method = "random"))
CDR3.overlap$overlap <- ifelse(CDR3.overlap$cdr3 %in% tcr.repertoire.data$CDR3, "yes", "no")

# stats
length(tcr.repertoire.data[which(tcr.repertoire.data$tcr_chain == "alpha"), "CDR3"])
length(which(tcr.repertoire.data[which(tcr.repertoire.data$tcr_chain == "alpha"), "CDR3"] %in%
  clonotypes$cdr3[which(clonotypes$chain == "TRA")]))
length(unique(clonotypes$cdr3[which(clonotypes$chain == "TRA")]))

length(tcr.repertoire.data[which(tcr.repertoire.data$tcr_chain == "beta"), "CDR3"])
length(which(tcr.repertoire.data[which(tcr.repertoire.data$tcr_chain == "beta"), "CDR3"] %in%
  clonotypes$cdr3[which(clonotypes$chain == "TRB")]))
length(unique(clonotypes$cdr3[which(clonotypes$chain == "TRB")]))

tcr.repertoire$overlap <- ifelse(tcr.repertoire$CDR3 %in% clonotypes$cdr3, "yes", "no")

boo <- tcr.repertoire %>%
  dplyr::select(sample, tcr_chain, overlap, frequency) %>%
  group_by(sample, tcr_chain, overlap) %>%
  summarize(freq = sum(frequency))
ggplot(boo[which(boo$overlap == "yes"), ], aes(x = sample, y = 100 * freq)) +
  geom_line(aes(group = tcr_chain, color = tcr_chain)) +
  scale_x_discrete(labels = c("Baseline", "C2D1", "C3D1", "C7D1", "Pause", "C2D1", "C3D1")) +
  scale_y_continuous("% overlap TCR reads PB - eAML", limits = c(0, 100)) +
  scale_color_manual(values = c("alpha" = "lightblue", "beta" = "darkblue")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230206_TCR_overlap_longitudinal.svg", width = 1.7, height = 2)

# p values
stats.df <- tcr.repertoire %>%
  group_by(sample, tcr_chain, overlap) %>%
  tally()

p.values.df <- data.frame()
for (chain in c("alpha", "beta")) {
  for (second.sample in c("44B", "44C", "44D", "44E", "44F", "44G")) {
    p.value <- fisher.test(matrix(c(
      stats.df$n[which(stats.df$sample == "44A" & stats.df$tcr_chain == "alpha" & stats.df$overlap == "yes")],
      stats.df$n[which(stats.df$sample == "44A" & stats.df$tcr_chain == "alpha" & stats.df$overlap == "no")],
      stats.df$n[which(stats.df$sample == second.sample & stats.df$tcr_chain == chain & stats.df$overlap == "yes")],
      stats.df$n[which(stats.df$sample == second.sample & stats.df$tcr_chain == chain & stats.df$overlap == "no")]
    ), ncol = 2))$p.value
    p.values.df <- rbind(p.values.df, data.frame(
      second.sample = second.sample,
      tcr.chain = chain,
      p.value = p.value
    ))
  }
}

p.values.df <- data.frame()
for (chain in c("alpha", "beta")) {
  for (sample.pair in list(
    c("44A", "44B"),
    c("44B", "44C"),
    c("44C", "44D"),
    c("44D", "44E"),
    c("44E", "44F"),
    c("44F", "44G")
  )) {
    p.value <- fisher.test(matrix(c(
      stats.df$n[which(stats.df$sample == sample.pair[1] & stats.df$tcr_chain == "alpha" & stats.df$overlap == "yes")],
      stats.df$n[which(stats.df$sample == sample.pair[1] & stats.df$tcr_chain == "alpha" & stats.df$overlap == "no")],
      stats.df$n[which(stats.df$sample == sample.pair[2] & stats.df$tcr_chain == chain & stats.df$overlap == "yes")],
      stats.df$n[which(stats.df$sample == sample.pair[2] & stats.df$tcr_chain == chain & stats.df$overlap == "no")]
    ), ncol = 2))$p.value
    p.values.df <- rbind(p.values.df, data.frame(
      first.sample = sample.pair[1],
      second.sample = sample.pair[2],
      tcr.chain = chain,
      p.value = p.value
    ))
  }
}



boo <- tcr.repertoire %>%
  group_by(sample, tcr_chain) %>%
  summarize(diversity = vegan::diversity(UMI) / log(length(UMI)))
ggplot(boo, aes(x = sample, y = diversity)) +
  geom_line(aes(group = tcr_chain, color = tcr_chain)) +
  scale_x_discrete(labels = c("Baseline", "C2D1", "C3D1", "C7D1", "Pause", "C2D1", "C3D1")) +
  scale_y_continuous("normalized Shannon index") +
  scale_color_manual(values = c("alpha" = "lightblue", "beta" = "darkblue")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230206_shannon_index_longitudinal.svg", width = 1.7, height = 2)


tcr.repertoire.filtered <- tcr.repertoire %>%
  group_by(tcr_chain, sample, UMI) %>%
  filter(n() > 1) %>%
  ungroup()

# remove a few insignificant duplicates
duplicates <- tcr.repertoire %>%
  dplyr::summarise(n = n(), .by = c(CDR3, sample)) %>%
  dplyr::filter(n > 1L)
tcr.repertoire <- tcr.repertoire[-which(tcr.repertoire$CDR3 %in% duplicates$CDR3), ]

tcr.repertoire.wide.alpha <- tidyr::pivot_wider(tcr.repertoire[which(tcr.repertoire$tcr_chain == "alpha"), ], id_cols = CDR3, names_from = sample, values_from = UMI) %>%
  as.data.frame()

tcr.repertoire.wide.beta <- tidyr::pivot_wider(tcr.repertoire[which(tcr.repertoire$tcr_chain == "beta"), ], id_cols = CDR3, names_from = sample, values_from = UMI) %>%
  as.data.frame()

# calculate p values
for (i in seq(2, 7)) {
  message(ecolTest::Hutcheson_t_test(
    as.numeric(tcr.repertoire.wide.alpha[, i]),
    as.numeric(tcr.repertoire.wide.alpha[, i + 1])
  )$p.value)
}
for (i in seq(2, 7)) {
  message(ecolTest::Hutcheson_t_test(
    as.numeric(tcr.repertoire.wide.beta[, i]),
    as.numeric(tcr.repertoire.wide.beta[, i + 1])
  )$p.value)
}


clonotypes$sample <- stringr::str_split_fixed(clonotypes$barcode, pattern = "_", n = 2)[, 1]
boo <- clonotypes %>%
  filter(cdr3 != "") %>%
  group_by(sample, chain) %>%
  summarize(
    cdr3.overlap = length(which(overlap == "yes")),
    cdr3.non.overlap = length(which(overlap == "no"))
  )
boo$ratio <- boo$cdr3.overlap / (boo$cdr3.overlap + boo$cdr3.non.overlap)

ggplot(
  CDR3.overlap[which(CDR3.overlap$chain == "TRA"), ],
  aes(x = rank, y = n)
) +
  geom_line(color = "grey") +
  geom_point(size = 1, aes(color = overlap)) +
  scale_x_log10() +
  scale_y_log10("cells") +
  scale_color_manual(values = c("yes" = "blue", "no" = "firebrick")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_scTCR_alpha.svg", width = 1.5, height = 1.5)

ggplot(
  CDR3.overlap[which(CDR3.overlap$chain == "TRB"), ],
  aes(x = rank, y = n)
) +
  geom_line(color = "grey") +
  geom_point(size = 1, aes(color = overlap)) +
  scale_x_log10() +
  scale_y_log10("cells") +
  scale_color_manual(values = c("yes" = "blue", "no" = "firebrick")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_scTCR_beta.svg", width = 1.5, height = 1.5)


tcr.repertoire.data <- tcr.repertoire.data %>%
  group_by(tcr_chain) %>%
  mutate(rank = rank(-UMI, ties.method = "random"))
tcr.repertoire.data$overlap <- ifelse(tcr.repertoire.data$CDR3 %in% clonotypes$cdr3, "yes", "no")

p <- ggplot(
  tcr.repertoire.data[which(tcr.repertoire.data$tcr_chain == "alpha"), ],
  aes(x = rank, y = UMI)
) +
  ggrastr::rasterize(geom_line(color = "grey"), dpi = 600) +
  ggrastr::rasterize(geom_point(size = 0.5, aes(color = overlap)), dpi = 600) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("yes" = "blue", "no" = "firebrick")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_bulk_TCR_alpha.svg", width = 1.5, height = 1.5, plot = p)

p <- ggplot(
  tcr.repertoire.data[which(tcr.repertoire.data$tcr_chain == "beta"), ],
  aes(x = rank, y = UMI)
) +
  ggrastr::rasterize(geom_line(color = "grey"), dpi = 600) +
  ggrastr::rasterize(geom_point(size = 0.5, aes(color = overlap)), dpi = 600) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("yes" = "blue", "no" = "firebrick")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_bulk_TCR_beta.svg", width = 1.5, height = 1.5, plot = p)

# donor vs. recipient
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.1", eAML.genotype$barcode) & eAML.genotype$assignment == "0")]) # recipient
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.1", eAML.genotype$barcode) & eAML.genotype$assignment == "1")]) # donor
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.2", eAML.genotype$barcode) & eAML.genotype$assignment == "0")]) # recipient
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.2", eAML.genotype$barcode) & eAML.genotype$assignment == "1")]) # donor
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.3", eAML.genotype$barcode) & eAML.genotype$assignment == "0")]) # recipient
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.3", eAML.genotype$barcode) & eAML.genotype$assignment == "1")]) # donor
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.4", eAML.genotype$barcode) & eAML.genotype$assignment == "0")]) # donor
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.4", eAML.genotype$barcode) & eAML.genotype$assignment == "1")]) # recipient
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.5", eAML.genotype$barcode) & eAML.genotype$assignment == "0")]) # recipient
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.5", eAML.genotype$barcode) & eAML.genotype$assignment == "1")]) # donor
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.6", eAML.genotype$barcode) & eAML.genotype$assignment == "0")]) # donor
DimPlot(eAML.data, cells.highlight = eAML.genotype$barcode[which(grepl("eAML1.6", eAML.genotype$barcode) & eAML.genotype$assignment == "1")]) # recipient

eAML.genotype$chimerism <- "none"
eAML.genotype$chimerism[which(grepl("eAML1.1", eAML.genotype$barcode) & eAML.genotype$assignment == "0")] <- "recipient"
eAML.genotype$chimerism[which(grepl("eAML1.1", eAML.genotype$barcode) & eAML.genotype$assignment == "1")] <- "donor"
eAML.genotype$chimerism[which(grepl("eAML1.2", eAML.genotype$barcode) & eAML.genotype$assignment == "0")] <- "recipient"
eAML.genotype$chimerism[which(grepl("eAML1.2", eAML.genotype$barcode) & eAML.genotype$assignment == "1")] <- "donor"
eAML.genotype$chimerism[which(grepl("eAML1.3", eAML.genotype$barcode) & eAML.genotype$assignment == "0")] <- "recipient"
eAML.genotype$chimerism[which(grepl("eAML1.3", eAML.genotype$barcode) & eAML.genotype$assignment == "1")] <- "donor"
eAML.genotype$chimerism[which(grepl("eAML1.4", eAML.genotype$barcode) & eAML.genotype$assignment == "0")] <- "donor"
eAML.genotype$chimerism[which(grepl("eAML1.4", eAML.genotype$barcode) & eAML.genotype$assignment == "1")] <- "recipient"
eAML.genotype$chimerism[which(grepl("eAML1.5", eAML.genotype$barcode) & eAML.genotype$assignment == "0")] <- "recipient"
eAML.genotype$chimerism[which(grepl("eAML1.5", eAML.genotype$barcode) & eAML.genotype$assignment == "1")] <- "donor"
eAML.genotype$chimerism[which(grepl("eAML1.6", eAML.genotype$barcode) & eAML.genotype$assignment == "0")] <- "donor"
eAML.genotype$chimerism[which(grepl("eAML1.6", eAML.genotype$barcode) & eAML.genotype$assignment == "1")] <- "recipient"
write.table(eAML.genotype, sep = "\t", quote = F, file = "./data/scRNA_seq/souporcell/20230119_eAML1_genotype.csv")

eAML.data$chimerism <- eAML.genotype[colnames(eAML.data), "chimerism"]
DimPlot(eAML.data, group.by = "chimerism", cols = c("donor" = "orange", "recipient" = "purple")) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_donor_recipient.png", width = 4, height = 4, dpi = 600)
FeaturePlot(eAML.data, "XIST", order = T) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank()) +
  scale_color_gradientn(colors = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_XIST.png", width = 4, height = 4, dpi = 600)
FeaturePlot(eAML.data, "RPS4Y1", order = T) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank()) +
  scale_color_gradientn(colors = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./scRNA_seq/figures/eAML1/20230119_eAML1_RPS4Y1.png", width = 4, height = 4, dpi = 600)

### T cells
eAML.Tcell <- subset(eAML.data, manual.cluster == "T cell")
eAML.Tcell <- NormalizeData(eAML.Tcell)
eAML.Tcell <- FindVariableFeatures(eAML.Tcell)
eAML.Tcell <- ScaleData(eAML.Tcell)
eAML.Tcell <- RunPCA(eAML.Tcell)
eAML.Tcell <- RunUMAP(eAML.Tcell, dims = 1:30)
eAML.Tcell <- FindNeighbors(eAML.Tcell)
eAML.Tcell <- FindClusters(eAML.Tcell, resolution = 0.3)

DimPlot(eAML.Tcell, cells.highlight = colnames(eAML.Tcell)[which(!is.na(eAML.Tcell$TCR_clonotype_id))])
DimPlot(eAML.Tcell, cells.highlight = clonotypes$barcode[which(clonotypes$cdr3 %in% tcr.repertoire.data$CDR3)], cols.highlight = "blue")
DimPlot(eAML.Tcell, cells.highlight = clonotypes$barcode[which(!clonotypes$cdr3 %in% tcr.repertoire.data$CDR3)], cols.highlight = "firebrick")

eAML.Tcell.markers <- FindAllMarkers(eAML.Tcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(eAML.Tcell.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC))

eAML.Tcell$manual.cluster <- "NK"
eAML.Tcell$manual.cluster[which(eAML.Tcell$seurat_clusters %in% c(0, 4, 5))] <- "CD8+ T cell"
eAML.Tcell$manual.cluster[which(eAML.Tcell$seurat_clusters %in% c(1, 2))] <- "CD4+ T cell"
eAML.Tcell$manual.cluster[which(eAML.Tcell$seurat_clusters %in% c(6))] <- "Treg"
eAML.Tcell$manual.cluster[which(eAML.Tcell$seurat_clusters %in% c(7, 8))] <- "CXCL13+ CD4+ T cell"
eAML.Tcell$manual.cluster <- factor(eAML.Tcell$manual.cluster, levels = names(eAML1.Tcell.colors))
DimPlot(eAML.Tcell, group.by = "manual.cluster", label = T)

eAML.Tcell <- ScaleData(eAML.Tcell, features = unique(c(Tcell.genes, "CXCL13", FindVariableFeatures(eAML.Tcell))))

saveRDS(eAML.Tcell, file = "./data/scRNA_seq/objects/20230120_eAML1_Tcells.rds")
eAML.Tcell <- readRDS(file = "./data/scRNA_seq/objects/20230120_eAML1_Tcells.rds")

VlnPlot(eAML.Tcell, "PDCD1") +
  scale_fill_manual(values = eAML1.Tcell.colors) +
  scale_y_continuous("PDCD1") +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    plot.title = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20241219_eAML1_Tcells_PDCD1.svg", width = 2, height = 3)

VlnPlot(eAML.Tcell, "CTLA4") +
  scale_fill_manual(values = eAML1.Tcell.colors) +
  scale_y_continuous("CTLA4") +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    plot.title = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20241219_eAML1_Tcells_CTLA4.svg", width = 2, height = 3)


eAML.Tcell$overlap <- "no.overlap"
eAML.Tcell$overlap[clonotypes$barcode[which(clonotypes$cdr3 %in% tcr.repertoire.data$CDR3)]] <- "overlap"
eAML.Tcell$overlap[which(eAML.Tcell$manual.cluster == "NK")] <- NA

ha <- columnAnnotation(
  celltype = eAML.Tcell$manual.cluster,
  overlap = eAML.Tcell$overlap,
  sample = eAML.Tcell$orig.ident,
  col = list(
    "celltype" = eAML1.Tcell.colors,
    "overlap" = c("overlap" = "blue", "no.overlap" = "firebrick"),
    "sample" = c(
      "eAML1.1" = RColorBrewer::brewer.pal(name = "Greens", n = 7)[2],
      "eAML1.2" = RColorBrewer::brewer.pal(name = "Greens", n = 7)[3],
      "eAML1.3" = RColorBrewer::brewer.pal(name = "Greens", n = 7)[4],
      "eAML1.4" = RColorBrewer::brewer.pal(name = "Greens", n = 7)[5],
      "eAML1.5" = RColorBrewer::brewer.pal(name = "Greens", n = 7)[6],
      "eAML1.6" = RColorBrewer::brewer.pal(name = "Greens", n = 7)[7]
    )
  ), na_col = "white",
  simple_anno_size = unit(5, "pt"), border = T, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 7, fontface = "bold")
)
df <- as.data.frame(GetAssayData(eAML.Tcell)[Tcell.genes, ])

col_fun <- circlize::colorRamp2(breaks = seq(0, 4, 4 / 8), colors = BuenColors::jdb_palette(name = "brewer_purple"))
svglite::svglite("./scRNA_seq/figures/eAML1/20230120_eAML1_Tcells_heatmap.svg", width = 7, height = 4.2)
Heatmap(df,
  column_split = eAML.Tcell$manual.cluster, top_annotation = ha, cluster_rows = F, cluster_columns = F,
  show_column_names = F, row_names_side = "left", row_names_gp = gpar(fontsize = 7, fontface = "italic"), column_title_gp = gpar(fontsize = 8), border = T,
  use_raster = T, raster_quality = 10, col = col_fun
)
dev.off()

Idents(eAML.Tcell) <- "manual.cluster"
df <- as.data.frame(prop.table(table(Idents(eAML.Tcell), eAML.Tcell$orig.ident), margin = 2))
df$Var1 <- factor(df$Var1, levels = names(eAML1.Tcell.colors))

ggplot(df, aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col() +
  scale_x_discrete(labels = c("pre-Nivo", "post-Nivo", "post-Ipi 1", "post-Ipi 2", "inflammation", "progression")) +
  scale_y_continuous("% cells") +
  scale_fill_manual("celltype", values = eAML1.Tcell.colors) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230120_eAML1_Tcell_celltypes_kinetics.svg", width = 3.5, height = 2)

# overlap between PB and eAML
df <- data.frame(
  barcode = colnames(eAML.Tcell),
  overlap = eAML.Tcell$overlap,
  manual.cluster = eAML.Tcell$manual.cluster
)

df <- as.data.frame(df %>% filter(manual.cluster != "NK") %>%
  group_by(manual.cluster, overlap) %>%
  summarize(count = n()) %>%
  tidyr::pivot_wider(names_from = overlap, values_from = count))
df$no.overlap.freq <- df$no.overlap / (df$no.overlap + df$overlap)

# p values

fisher.test(matrix(c(
  df$no.overlap[which(df$manual.cluster == "CD4+ T cell")],
  df$overlap[which(df$manual.cluster == "CD4+ T cell")],
  df$no.overlap[which(df$manual.cluster == "CXCL13+ CD4+ T cell")],
  df$overlap[which(df$manual.cluster == "CXCL13+ CD4+ T cell")]
), ncol = 2))

fisher.test(matrix(c(
  df$no.overlap[which(df$manual.cluster == "CD8+ T cell")],
  df$overlap[which(df$manual.cluster == "CD8+ T cell")],
  df$no.overlap[which(df$manual.cluster == "CXCL13+ CD4+ T cell")],
  df$overlap[which(df$manual.cluster == "CXCL13+ CD4+ T cell")]
), ncol = 2))

fisher.test(matrix(c(
  df$no.overlap[which(df$manual.cluster == "Treg")],
  df$overlap[which(df$manual.cluster == "Treg")],
  df$no.overlap[which(df$manual.cluster == "CXCL13+ CD4+ T cell")],
  df$overlap[which(df$manual.cluster == "CXCL13+ CD4+ T cell")]
), ncol = 2))

ggplot(df, aes(x = manual.cluster, y = 100 * no.overlap.freq, fill = manual.cluster)) +
  geom_col() +
  scale_y_continuous("% non-overlapping T cells ") +
  scale_fill_manual("celltype", values = eAML1.Tcell.colors) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230120_eAML1_Tcells_overlap.svg", width = 3, height = 3)

### combined clonotypes
source("./scRNA_seq/R/00_read_clonotypes.R")

eAML.Tcell <- read_clonotypes("./data/scRNA_seq/20230207_TCR_eAML1.xlsx", AML.combined = eAML.Tcell, TCR.clones.file = "./data/scRNA_seq/20230207_TCR_eAML1.csv")

### visualize TCR repertoire
vdjdb <- data.table::fread("/Users/shaka87/dfci/vdjdb/vdjdb_full.txt") %>% filter(vdjdb.score > 1)
tcr.repertoire.data$vdjdb <- ifelse(tcr.repertoire.data$CDR3 %in% vdjdb$cdr3.alpha | tcr.repertoire.data$CDR3 %in% vdjdb$cdr3.beta, "yes", "no")
tcr.repertoire.data$overlap <- ifelse(tcr.repertoire.data$CDR3 %in% clonotypes$cdr3, "yes", "no")
clonotypes$vdjdb <- ifelse(clonotypes$cdr3 %in% vdjdb$cdr3.alpha | clonotypes$cdr3 %in% vdjdb$cdr3.beta, "yes", "no")

eAML.Tcell$vdjdb <- "no"
eAML.Tcell$vdjdb[clonotypes$barcode[which(clonotypes$vdjdb == "yes")]] <- "yes"
eAML.Tcell$overlap <- "no"
eAML.Tcell$overlap[clonotypes$barcode[which(clonotypes$overlap == "yes")]] <- "yes"

df <- data.frame(
  bc = colnames(eAML.Tcell),
  sample = eAML.Tcell$orig.ident,
  manual.cluster = eAML.Tcell$manual.cluster,
  overlap = eAML.Tcell$overlap
)

df.stats <- df %>%
  group_by(sample, manual.cluster) %>%
  summarize(
    cells.overlap = length(which(overlap == "yes")),
    cells.non.overlap = length(which(overlap == "no"))
  )
df.stats$cells <- df.stats$cells.overlap + df.stats$cells.non.overlap
df.stats$overlap.freq <- df.stats$cells.overlap / df.stats$cells

# overlap PB and eAML
ggplot(df.stats[which(df.stats$manual.cluster != "NK"), ], aes(x = cells, y = 100 * overlap.freq, color = manual.cluster)) +
  geom_point(size = 1) +
  scale_x_continuous("cells") +
  scale_y_continuous("%detected in PB", limits = c(0, 100)) +
  scale_color_manual(values = eAML1.Tcell.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/eAML1/20230207_TCR_overlap.svg", width = 2, height = 2)


# TCR repertoire as barplots
df <- data.frame(
  bc = colnames(eAML.Tcell),
  sample = eAML.Tcell$orig.ident,
  manual.cluster = eAML.Tcell$manual.cluster,
  CloneID = eAML.Tcell$CloneID
)
df <- df[-which(df$manual.cluster == "NK"), ]


boo <- df %>%
  filter(sample == "eAML1.1") %>%
  group_by(CloneID) %>%
  summarize(
    CD4 = length(which(manual.cluster == "CD4+ T cell")),
    CD4.ex = length(which(manual.cluster == "CXCL13+ CD4+ T cell")),
    CD8 = length(which(manual.cluster == "CD8+ T cell")),
    Treg = length(which(manual.cluster == "Treg")),
    cells = length(CloneID)
  )
boo <- boo[-which(is.na(boo$CloneID)), ]
boo <- boo[order(-boo$cells), ]
boo$rank <- seq(1, nrow(boo))
boo <- boo %>% tidyr::pivot_longer(cols = c("CD8", "CD4", "CD4.ex", "Treg"), names_to = "manual.cluster", values_to = "value")
boo <- boo[-which(boo$value == 0), ]
boo <- boo[c(
  which(boo$cells > 1),
  which(boo$cells == 1 & boo$manual.cluster == "CD4"),
  which(boo$cells == 1 & boo$manual.cluster == "CD4.ex"),
  which(boo$cells == 1 & boo$manual.cluster == "Treg"),
  which(boo$cells == 1 & boo$manual.cluster == "CD8")
), ]
boo$rank[which(boo$cells == 1)] <- seq(min(boo$rank[which(boo$cells == 1)]), max(boo$rank[which(boo$cells == 1)]))

ggplot(boo, aes(x = rank, y = value)) +
  geom_col(aes(fill = manual.cluster)) +
  scale_y_continuous("clonal T cell expansion") +
  scale_fill_manual(values = eAML1.Tcell.colors2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230207_TCR_repertoire_eAML1.1.svg", width = 2, height = 1.5)

boo <- df %>%
  filter(sample == "eAML1.2") %>%
  group_by(CloneID) %>%
  summarize(
    CD4 = length(which(manual.cluster == "CD4+ T cell")),
    CD4.ex = length(which(manual.cluster == "CXCL13+ CD4+ T cell")),
    CD8 = length(which(manual.cluster == "CD8+ T cell")),
    Treg = length(which(manual.cluster == "Treg")),
    cells = length(CloneID)
  )
boo <- boo[-which(is.na(boo$CloneID)), ]
boo <- boo[order(-boo$cells), ]
boo$rank <- seq(1, nrow(boo))
boo <- boo %>% tidyr::pivot_longer(cols = c("CD8", "CD4", "CD4.ex", "Treg"), names_to = "manual.cluster", values_to = "value")
boo <- boo[-which(boo$value == 0), ]
boo <- boo[c(
  which(boo$cells > 1),
  which(boo$cells == 1 & boo$manual.cluster == "CD4"),
  which(boo$cells == 1 & boo$manual.cluster == "CD4.ex"),
  which(boo$cells == 1 & boo$manual.cluster == "Treg"),
  which(boo$cells == 1 & boo$manual.cluster == "CD8")
), ]
boo$rank[which(boo$cells == 1)] <- seq(min(boo$rank[which(boo$cells == 1)]), max(boo$rank[which(boo$cells == 1)]))

ggplot(boo, aes(x = rank, y = value)) +
  geom_col(aes(fill = manual.cluster)) +
  scale_y_continuous("clonal T cell expansion") +
  scale_fill_manual(values = eAML1.Tcell.colors2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230207_TCR_repertoire_eAML1.2.svg", width = 2, height = 1.5)


boo <- df %>%
  filter(sample == "eAML1.5") %>%
  group_by(CloneID) %>%
  summarize(
    CD4 = length(which(manual.cluster == "CD4+ T cell")),
    CD4.ex = length(which(manual.cluster == "CXCL13+ CD4+ T cell")),
    CD8 = length(which(manual.cluster == "CD8+ T cell")),
    Treg = length(which(manual.cluster == "Treg")),
    cells = length(CloneID)
  )
boo <- boo[-which(is.na(boo$CloneID)), ]
boo <- boo[order(-boo$cells), ]
boo$rank <- seq(1, nrow(boo))
boo <- boo %>% tidyr::pivot_longer(cols = c("CD8", "CD4", "CD4.ex", "Treg"), names_to = "manual.cluster", values_to = "value")
boo <- boo[-which(boo$value == 0), ]
boo <- boo[c(
  which(boo$cells > 1),
  which(boo$cells == 1 & boo$manual.cluster == "CD4"),
  which(boo$cells == 1 & boo$manual.cluster == "CD4.ex"),
  which(boo$cells == 1 & boo$manual.cluster == "Treg"),
  which(boo$cells == 1 & boo$manual.cluster == "CD8")
), ]
boo$rank[which(boo$cells == 1)] <- seq(min(boo$rank[which(boo$cells == 1)]), max(boo$rank[which(boo$cells == 1)]))

ggplot(boo, aes(x = rank, y = value)) +
  geom_col(aes(fill = manual.cluster)) +
  scale_y_continuous("clonal T cell expansion") +
  scale_fill_manual(values = eAML1.Tcell.colors2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave("./scRNA_seq/figures/eAML1/20230207_TCR_repertoire_eAML1.5.svg", width = 2, height = 1.5)

# overview of TCR repertoire as heatmap
df <- data.frame(
  bc = colnames(eAML.Tcell),
  sample = eAML.Tcell$orig.ident,
  manual.cluster = eAML.Tcell$manual.cluster,
  CloneID = eAML.Tcell$CloneID
)
df <- df[-which(is.na(df$CloneID)), ]
df <- df[-which(df$manual.cluster == "NK"), ]

df.stats <- df %>%
  group_by(sample, CloneID) %>%
  summarize(n = n())
df.stats <- df.stats %>%
  group_by(sample) %>%
  mutate(freq = n / length(sample)) %>%
  dplyr::select(-n)
df.stats <- df.stats %>% tidyr::pivot_wider(names_from = "sample", values_from = "freq")
df.stats <- df.stats[order(-df.stats$eAML1.1), ]
col_fun <- circlize::colorRamp2(breaks = seq(-3, -1, 2 / 9), colors = c("white", BuenColors::jdb_palette(name = "solar_rojos")))

svglite::svglite("./scRNA_seq/figures/eAML1/20230207_TCR_repertoire_heatmap.svg", width = 2, height = 2)
Heatmap(log10(df.stats[, -1]),
  show_row_names = F, show_column_names = T, cluster_rows = F,
  column_labels = c("pre-Nivo", "post-Nivo", "post-Ipi 1", "post-Ipi 2", "response"),
  column_names_gp = gpar(fontsize = 10),
  cluster_columns = F, na_col = "white", border = T, col = col_fun
)
dev.off()
