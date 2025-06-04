# data on macrophages 

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(Seurat)

source('./scRNA_seq/R/00_eAML.colors.R')

eAML1 = readRDS('./data/scRNA_seq/objects/20230117_eAML1.rds')

# https://pmc.ncbi.nlm.nih.gov/articles/PMC4045180/

genes = c('CD14', 'FCGR1A',  'CX3CR1', 'CSF1R', 'ITGAM', 'ITGAX','SIGLEC1', 'CLEC7A', 'CLEC10A', 'CD80', 'CD86', 'HLA-DRA', 'CD68', 'CD40', 'CD163', 'MRC1', 'CD209')

eAML1.macrophages = subset(eAML1, manual.cluster == 'macrophage')
eAML1.macrophages = RunPCA(eAML1.macrophages)
eAML1.macrophages = FindNeighbors(eAML1.macrophages, dims = 1:7)
eAML1.macrophages = FindClusters(eAML1.macrophages)
eAML1.macrophages = RunUMAP(eAML1.macrophages, dims = 1:7)
eAML1.macrophages = ScaleData(eAML1.macrophages, features = genes)

mat = GetAssayData(eAML1.macrophages, layer = 'scale.data')[genes,] %>% as.matrix()
col_fun = circlize::colorRamp2(breaks = seq(-2, 2.5, 4.5/8), BuenColors::jdb_palette(name = 'brewer_celsius', n=9))
ha = columnAnnotation(sample = eAML1.macrophages$orig.ident)
svglite::svglite('./scRNA_seq/figures/combined/heatmap/20250502_macrophage_phenotype.svg', width = 4.5, height = 3)
Heatmap(mat, show_column_names = F, col = col_fun, 
        cluster_columns = F, border = T, row_names_side = 'left', row_names_gp = gpar(fontsize=10, fontface='italic'),
        cluster_rows = F, top_annotation = ha, column_split = eAML1.macrophages$orig.ident)
dev.off()
