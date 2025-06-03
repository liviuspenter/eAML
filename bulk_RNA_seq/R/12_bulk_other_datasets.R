library(biomaRt)
library(ComplexHeatmap)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(stringr)
library(tximport)
library(GenomicFeatures)

# mapping of transcripts to genes
txdb <- makeTxDbFromGFF("./data/GRCh38/gencode.v43.annotation.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- ensembldb::select(txdb, k, "GENEID", "TXNAME")

tx2gene$TXNAME <- gsub(tx2gene$TXNAME, pattern = "\\..*", replacement = "")
tx2gene$GENEID <- gsub(tx2gene$GENEID, pattern = "\\..*", replacement = "")

# output from salmon
files <- list.files(path = "./data/bulk/quants/", pattern = "SRR|CIMAC")
files_import <- paste0("./data/bulk/quants/", files, "/quant.sf")

files <- files[which(file.exists(files_import))]
files_import <- files_import[which(file.exists(files_import))]

# read salmon output
mat_gse <- tximport(files_import,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE, ignoreTxVersion = T
)

TPM.mat <- as.data.frame(mat_gse$abundance)
colnames(TPM.mat) <- files

# get gene names
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(TPM.mat)
genes <- stringr::str_split_fixed(genes, pattern = "\\.", n = 2)[, 1]
gene_IDs <- getBM(
  filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = genes, mart = mart
)
gene_IDs <- as.data.frame(gene_IDs %>% distinct(ensembl_gene_id, .keep_all = T))
rownames(gene_IDs) <- gene_IDs$ensembl_gene_id

TPM.mat$ensembl_gene_id <- rownames(TPM.mat)
TPM.mat$gene <- gene_IDs[genes, "hgnc_symbol"]
TPM.mat <- TPM.mat[, c("ensembl_gene_id", "gene", files)]

write.csv2(TPM.mat, file = "./data/bulk/20250416_TPM.csv", quote = F)
write.table(TPM.mat[, -1], file = "./data/bulk/20250416_TPM_cibertsort_input.csv", quote = F, row.names = F, dec = ".", sep = "\t")

TPM.mat <- data.table::fread(file = "./data/bulk/20250416_TPM.csv", header = T, dec = ",") %>% as.data.frame()

metadata <- as.data.frame(readxl::read_excel("./data/bulk/20250416_eAML_bulk_metadata.xlsx", sheet = 1))
rownames(metadata) <- metadata$library

TPM.mat <- TPM.mat[, c("ensembl_gene_id", "gene", intersect(metadata$library, colnames(TPM.mat)))]

genes <- c(
  "CD3D", "CD4", "CD8A", "NCAM1", "SELL", "CCR7", "IL7R", "TCF7", "CD28", "FAS", "CD27", "ITGAE",
  "CTLA4", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "ICOS", "ENTPD1", "BTLA", "KLRG1", "KLRB1", "CD38", "CD69", "FOXP3",
  "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "NKG7", "GNLY"
)

mat <- TPM.mat[which(TPM.mat$gene %in% genes), metadata$library[which(metadata$status %in% c("skin", "eAML", "GvHD", "LS"))]]
# rownames(mat) = gene_IDs[as.integer(rownames(mat)), 'hgnc_symbol']
rownames(mat) <- TPM.mat$gene[which(TPM.mat$gene %in% genes)]

cols <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 9), colors = c("white", BuenColors::jdb_palette(name = "brewer_green")))
svglite::svglite("./bulk_RNA_seq/figures/combined/heatmap/20250416_eAML_skin_GvHD_MPR_LS_bulk_Tcell_genes.svg", width = 3, height = 2.9)
Heatmap(mat[genes, ],
  cluster_rows = F, cluster_columns = F, row_names_side = "left", row_names_gp = gpar(fontsize = 7, fontface = "italic"),
  show_column_names = F, col = cols,
  column_split = factor(metadata$status[which(metadata$status %in% c("skin", "eAML", "GvHD", "LS"))], levels = c("skin", "eAML", "GvHD", "LS")),
  border = T
)
dev.off()


### explore differences between scleroderma and eAML

# mapping of transcripts to genes
txdb <- makeTxDbFromGFF("./data/GRCh38/gencode.v43.annotation.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- ensembldb::select(txdb, k, "GENEID", "TXNAME")

samples <- list.files(path = "./data/bulk/quants/", full.names = T, pattern = "SRR|CIMAC")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./data/bulk/quants//", "")
files <- files[intersect(names(files), metadata$library[which(metadata$status %in% c("eAML", "LS"))])]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = F)

# get gene names
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(txi$counts)
genes <- stringr::str_split_fixed(genes, pattern = "\\.", n = 2)[, 1]
gene_IDs <- getBM(
  filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = genes, mart = mart
)
gene_IDs <- as.data.frame(gene_IDs %>% distinct(ensembl_gene_id, .keep_all = T))
rownames(gene_IDs) <- gene_IDs$ensembl_gene_id

sampletype <- factor(metadata[colnames(txi$counts), "status"])
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

#
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~sampletype)
dds <- DESeq(dds)

normalized_counts <- counts(dds, normalized = T)

# rld = rlog(dds, blind = T)

rld <- vst(dds, blind = T)

plotPCA(rld, intgroup = "sampletype")

contrast_oe <- c("sampletype", "LS", "eAML")
res_tableOE <- results(dds, contrast = contrast_oe, alpha = 0.05)
res_tableOE_shrunk <- lfcShrink(dds, coef = "sampletype_LS_vs_eAML", type = "ashr")

res_tableOE_tb <- res_tableOE_shrunk %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  as_tibble()

sigOE <- res_tableOE_tb %>% filter(padj < 0.05 & abs(log2FoldChange) > 2)
sigOE$gene.id <- gene_IDs[gsub(sigOE$gene, pattern = "\\..*", replacement = ""), "hgnc_symbol"]

res_tableOE_tb$gene.id <- gene_IDs[gsub(res_tableOE_tb$gene, pattern = "\\..*", replacement = ""), "hgnc_symbol"]

highlight.genes <- c("SELL", "CD38", "GZMK", "WT1", "RNASE2", "MPO", "CD33")

p <- ggplot() +
  ggrastr::rasterize(geom_point(
    data = sigOE, aes(x = -log2FoldChange, y = -log10(padj)),
    color = "grey", size = 0.5
  ), dpi = 600) +
  ggrastr::rasterize(geom_point(
    data = sigOE[which(sigOE$gene.id %in% highlight.genes), ], aes(x = -log2FoldChange, y = -log10(padj)),
    color = "firebrick", size = 0.5
  ), dpi = 600) +
  geom_label_repel(
    data = sigOE[which(sigOE$gene.id %in% highlight.genes), ], aes(x = -log2FoldChange, y = -log10(padj), label = gene.id),
    color = "black", label.size = 0, max.overlaps = 15, size = 3
  ) +
  theme_classic() +
  scale_x_continuous("log2FC", limits = c(-15, 15)) +
  scale_y_continuous("-log10(FDR)", limits = c(0, 30)) +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./bulk_RNA_seq/figures/combined/volcano/20250416_bulk_RNA_eAML_LS.svg", width = 3, height = 2.5, plot = p)


### CIBERSORTx
cibertsortx.data <- as.data.frame(data.table::fread("./data/bulk/20250416_CIBERSORTx_lm22.csv"))
rownames(cibertsortx.data) <- cibertsortx.data$Mixture
cibertsortx.data <- cibertsortx.data[metadata$library[which(metadata$status %in% c("eAML", "skin", "GvHD", "LS"))], ]
cibertsortx.data <- cibertsortx.data[, -1]
cibertsortx.data <- cibertsortx.data[, ncol(cibertsortx.data) - 4:ncol(cibertsortx.data)]

cibertsortx.data <- as.data.frame(t(cibertsortx.data))

ha <- columnAnnotation(
  tissue = metadata[colnames(cibertsortx.data), "status"],
  col = list("tissue" = c("eAML" = "orange", skin = "blue", "GvHD" = "red", "LS" = "brown")),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_side = "left",
  simple_anno_size = unit(5, "pt"), border = T
)

cols <- circlize::colorRamp2(breaks = seq(0, 8, 1), colors = BuenColors::jdb_palette(name = "samba_color"))
svglite::svglite("./bulk_RNA_seq/figures/combined/heatmap/20250416_eAML_skin_bulk.svg", width = 4.5, height = 3)
Heatmap(cibertsortx.data,
  cluster_columns = F, cluster_rows = T, col = cols,
  column_split = factor(metadata[colnames(cibertsortx.data), "status"], levels = c("skin", "eAML", "GvHD", "LS")), top_annotation = ha,
  row_names_side = "left", row_names_gp = gpar(fontsize = 8), show_column_names = F, border = T,
  column_title_gp = gpar(fontsize = 10)
) # , right_annotation = ha2)
dev.off()

# perform differential abundance analysis eAML vs LS
LS.samples <- which(colnames(cibertsortx.data) %in% metadata$library[which(metadata$status == "LS")])
eAML.samples <- which(colnames(cibertsortx.data) %in% metadata$library[which(metadata$status == "eAML")])

p.values <- apply(cibertsortx.data, MARGIN = 1, FUN = function(x) {
  t.test(x[LS.samples], x[eAML.samples])$p.value
})
cibertsortx.data <- cibertsortx.data[names(sort(p.values)), ]

df <- data.frame(
  log2fc = log2(rowMeans(cibertsortx.data[, LS.samples]) /
    rowMeans(cibertsortx.data[, eAML.samples])),
  p.values = p.values
)

df$significant <- "none"
df[which(df$log2fc > 1 & df$p.values < 0.05), "significant"] <- "LS"
df[which(df$log2fc < 1 & df$p.values < 0.05), "significant"] <- "eAML"
ggplot(df, aes(x = log2fc, y = -log10(p.values))) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = c(-1, 1)) +
  geom_point(aes(color = significant), size = 0.5) +
  geom_label_repel(
    data = df[which(abs(df$log2fc) > 1 & df$p.values < 0.05), ],
    aes(
      x = log2fc, y = -log10(p.values),
      label = rownames(df[which(abs(df$log2fc) > 1 & df$p.values < 0.05), ])
    ), label.size = 0, size = 3
  ) +
  scale_color_manual(values = c("eAML" = "orange", "LS" = "brown", "GvHD" = "red", "none" = "grey")) +
  scale_x_continuous("log2FC") +
  scale_y_continuous("-log10(p-value)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./bulk_RNA_seq/figures/combined/volcano/20250416_volcano_LS_eAML.svg", width = 2.5, height = 2.5)

# perform differential abundance analysis eAML vs GvHD plus LS
other.samples <- which(colnames(cibertsortx.data) %in% metadata$library[which(metadata$status %in% c("GvHD", "LS"))])
eAML.samples <- which(colnames(cibertsortx.data) %in% metadata$library[which(metadata$status == "eAML")])

p.values <- apply(cibertsortx.data, MARGIN = 1, FUN = function(x) {
  t.test(x[other.samples], x[eAML.samples])$p.value
})
cibertsortx.data <- cibertsortx.data[names(sort(p.values)), ]

df <- data.frame(
  log2fc = log2(rowMeans(cibertsortx.data[, other.samples]) /
    rowMeans(cibertsortx.data[, eAML.samples])),
  p.values = p.values
)

df$significant <- "none"
df[which(df$log2fc > 1 & df$p.values < 0.05), "significant"] <- "LS_GvHD"
df[which(df$log2fc < 1 & df$p.values < 0.05), "significant"] <- "eAML"
ggplot(df, aes(x = log2fc, y = -log10(p.values))) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = c(-1, 1)) +
  geom_point(aes(color = significant), size = 0.5) +
  geom_label_repel(
    data = df[which(abs(df$log2fc) > 1 & df$p.values < 0.05), ],
    aes(
      x = log2fc, y = -log10(p.values),
      label = rownames(df[which(abs(df$log2fc) > 1 & df$p.values < 0.05), ])
    ), label.size = 0, size = 3
  ) +
  scale_color_manual(values = c("eAML" = "orange", "LS" = "brown", "GvHD" = "red", "none" = "grey", "LS_GvHD" = "black")) +
  scale_x_continuous("log2FC") +
  scale_y_continuous("-log10(p-value)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./bulk_RNA_seq/figures/combined/volcano/20250416_volcano_LS_GvHD_eAML.svg", width = 2.5, height = 2.5)


df <- t(cibertsortx.data) %>% as.data.frame()
df$tissue <- factor(metadata[colnames(cibertsortx.data), "status"],
  levels = c("skin", "eAML", "GvHD", "LS")
)

ggplot(df, aes(x = tissue, y = `Macrophages M2`)) +
  geom_boxplot(aes(fill = tissue), outliers = F, color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(
    c("skin", "eAML"), c("skin", "LS"),
    c("eAML", "LS")
  ), step_increase = 0.15, textsize = 3, test = "t.test") +
  theme_classic() +
  scale_fill_manual(values = c("skin" = "blue", "eAML" = "orange", "LS" = "brown", "GvHD" = "red")) +
  scale_x_discrete(labels = c("healthy skin", "leukemia cutis", "GvHD", "scleroderma")) +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./bulk_RNA_seq/figures/combined/plots/20250416_macrophages_M2.svg", width = 1.5, height = 3)


ggplot(df, aes(x = tissue, y = `NK cells activated`)) +
  geom_boxplot(aes(fill = tissue), outliers = F) +
  geom_jitter(width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(
    c("skin", "eAML"), c("skin", "GvHD"), c("skin", "LS"),
    c("eAML", "LS")
  ), step_increase = 0.15, textsize = 3) +
  theme_classic() +
  scale_fill_manual(values = c("skin" = "blue", "eAML" = "orange", "LS" = "brown", "GvHD" = "red")) +
  scale_x_discrete(labels = c("skin", "eAML", "GvHD", "scleroderma")) +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./bulk_RNA_seq/figures/combined/plots/20250416_NK_activated.svg", width = 2, height = 3)

ggplot(df, aes(x = tissue, y = `NK cells resting`)) +
  geom_boxplot(aes(fill = tissue), outliers = F) +
  geom_jitter(width = 0.1, size = 0.5) +
  geom_signif(comparisons = list(
    c("skin", "eAML"), c("skin", "GvHD"), c("skin", "LS"),
    c("eAML", "LS")
  ), step_increase = 0.15, textsize = 3) +
  theme_classic() +
  scale_fill_manual(values = c("skin" = "blue", "eAML" = "orange", "LS" = "brown", "GvHD" = "red")) +
  scale_x_discrete(labels = c("skin", "eAML", "GvHD", "scleroderma")) +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./bulk_RNA_seq/figures/combined/plots/20250416_NK_resting.svg", width = 2, height = 3)
