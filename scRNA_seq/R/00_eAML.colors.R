eAML.combined.colors <- c(
  "myeloid" = BuenColors::jdb_palette(name = "corona")[2],
  "macrophage" = BuenColors::jdb_palette(name = "corona")[4],
  "T" = BuenColors::jdb_palette(name = "corona")[1],
  "NK" = BuenColors::jdb_palette(name = "corona")[10],
  "B cell" = BuenColors::jdb_palette(name = "corona")[3],
  "fibroblast" = BuenColors::jdb_palette(name = "corona")[6],
  "keratinocyte" = BuenColors::jdb_palette(name = "corona")[13],
  "endothelium" = BuenColors::jdb_palette(name = "corona")[5],
  "melanocyte" = BuenColors::jdb_palette(name = "corona")[11],
  "smooth muscle" = BuenColors::jdb_palette(name = "corona")[12],
  "erythroid" = BuenColors::jdb_palette(name = "corona")[8],
  "platelet" = BuenColors::jdb_palette(name = "corona")[15]
)


eAML.Tcell.colors <- c(
  "naive_CM" = "#06D6A0", "CD4+ T cell" = "#FFD116", "Treg" = "#F78C6B", "CD8+ T cell" = "#EF476F",
  "CD16+ NK" = "#118AB2", "NK" = "#073B4C", "doublet" = "grey90"
)

donor.recipient.colors <- c("recipient" = "purple", "donor" = "orange")

eAML1.colors <- c(
  "AML" = BuenColors::jdb_palette(name = "corona")[1],
  "T cell" = BuenColors::jdb_palette(name = "corona")[3],
  "B cell" = BuenColors::jdb_palette(name = "corona")[2],
  "keratinocyte" = BuenColors::jdb_palette(name = "corona")[7],
  "endothelium" = BuenColors::jdb_palette(name = "corona")[5],
  "macrophage" = BuenColors::jdb_palette(name = "corona")[4],
  "fibroblast" = BuenColors::jdb_palette(name = "corona")[6]
)

eAML1.colors.chimerism <- rep(c("orange", "purple"), 7)
names(eAML1.colors.chimerism) <- paste0(rep(names(eAML1.colors), each = 2), c(".donor", ".recipient"))

eAML1.Tcell.colors <- c(
  "CD8+ T cell" = "orange",
  "CD4+ T cell" = RColorBrewer::brewer.pal(name = "Blues", n = 4)[2],
  "CXCL13+ CD4+ T cell" = RColorBrewer::brewer.pal(name = "Blues", n = 4)[3],
  "Treg" = RColorBrewer::brewer.pal(name = "Blues", n = 4)[4],
  "NK" = "black"
)
eAML1.Tcell.colors2 <- eAML1.Tcell.colors
names(eAML1.Tcell.colors2) <- c("CD8", "CD4", "CD4.ex", "Treg", "NK")

eAML2.colors <- c(
  "AML" = BuenColors::jdb_palette(name = "corona")[1],
  "B cell" = BuenColors::jdb_palette(name = "corona")[2],
  "T cell" = BuenColors::jdb_palette(name = "corona")[3],
  "endothelium" = BuenColors::jdb_palette(name = "corona")[5],
  "fibroblast" = BuenColors::jdb_palette(name = "corona")[6],
  "keratinocyte" = BuenColors::jdb_palette(name = "corona")[7],
  "smooth muscle" = BuenColors::jdb_palette(name = "corona")[8]
)

eAML3.colors <- c(
  "AML" = BuenColors::jdb_palette(name = "corona")[1],
  "T cell" = BuenColors::jdb_palette(name = "corona")[3],
  "endothelium" = BuenColors::jdb_palette(name = "corona")[5],
  "fibroblast" = BuenColors::jdb_palette(name = "corona")[6],
  "keratinocyte" = BuenColors::jdb_palette(name = "corona")[7],
  "smooth muscle" = BuenColors::jdb_palette(name = "corona")[8],
  "melanocyte" = BuenColors::jdb_palette(name = "corona")[9],
  "gland" = BuenColors::jdb_palette(name = "corona")[10],
  "platelet" = BuenColors::jdb_palette(name = "corona")[11],
  "erythrocyte" = BuenColors::jdb_palette(name = "corona")[12]
)

Tcell.genes <- c(
  "CD3D", "CD3E", "CD3G", "CD8A", "CD4", "SELL", "CCR7", "IL7R", "CD27", "CD28", "FAS", "NCAM1", "B3GAT1",
  "PDCD1", "CTLA4", "TIGIT", "HAVCR2", "LAG3", "BTLA", "CD244", "KLRG1", "KLRB1",
  "ENTPD1", "ITGAE", "CD69", "IL2RA", "ICOS", "CXCL13",
  "TNFRSF4", "TNFRSF9", "HLA-DRA", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "IFNG", "FASLG", "TNF", "IL10",
  "TCF7", "EOMES", "TBX21", "TOX", "FOXP3"
)

tissue.colors <- c("BM" = "firebrick", "PB" = "red", "eAML" = "orange", "skin" = "blue")

BM.samples <- c("1007", "1010", "1012", "1016", "1019", "1022", "1026")
eAML.samples <- c("eAML1.1", "eAML1.2", "eAML1.3", "eAML1.4", "eAML1.5", "eAML2", "eAML3.1", "eAML3.2", "eAML4.1", "eAML4.2", "eAML5.1", "eAML5.2", "eAML6", "eAML7.1", "eAML7.2")
eAML.samples2 <- c("eAML1.1", "eAML1.2", "eAML1.3", "eAML1.4", "eAML1.5", "eAML1.6", "eAML2", "eAML3.1", "eAML3.2", "eAML4.1", "eAML4.2", "eAML7.1", "eAML7.2", "eAML5.1", "eAML5.2", "eAML6")
PB.samples <- c("eAML2PB", "eAML3PB", "eAML4PB.1", "eAML4PB.2", "eAML4PB.3", "eAML7PB")
skin.samples <- c("P112", "P115", "P116", "P121", "SKN8104902", "SKN8105200")

BM.samples.colors <- c("firebrick", "white", "firebrick", "white", "firebrick", "white", "firebrick")
names(BM.samples.colors) <- as.character(BM.samples)
eAML.samples.colors <- c("orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange")
names(eAML.samples.colors) <- eAML.samples
eAML.samples.colors2 <- c("orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange", "white", "orange", "white")
names(eAML.samples.colors2) <- eAML.samples2
PB.samples.colors <- c("red", "white", "red", "white", "red", "white")
names(PB.samples.colors) <- PB.samples
skin.samples.colors <- c("blue", "white", "blue", "white", "blue", "white")
names(skin.samples.colors) <- skin.samples
