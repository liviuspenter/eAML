library(ggplot2)

metadata <- as.data.frame(readxl::read_excel("./data/bulk/20230302_eAML_bulk_metadata.xlsx", sheet = 1))
rownames(metadata) <- metadata$library

df <- data.frame()
for (f in list.files(path = "./data/bulk/TRUST4/")) {
  library.name <- stringr::str_split_fixed(f, pattern = "_", n = 3)[, 1]

  boo <- data.table::fread(paste0("./data/bulk/TRUST4/", f))

  df <- rbind(df, data.frame(
    library = library.name,
    TRA = length(which(boo$V6 == "TRAC")),
    TRB = length(which(boo$V6 %in% c("TRBC1", "TRBC2"))),
    TRD = length(which(boo$V6 %in% c("TRDC"))),
    TRG = length(which(boo$V6 %in% c("TRGC2"))),
    IGH = length(which(grepl("IGH", boo$V6))),
    IGK = length(which(grepl("IGK", boo$V6))),
    IGL = length(which(grepl("IGL", boo$V6)))
  ))
}

df$sampletype <- metadata[df$library, "status"]

boo <- df %>%
  filter(!is.na(sampletype)) %>%
  tidyr::pivot_longer(names_to = "receptor", cols = c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"))
boo$receptor.sampletype <- factor(paste0(boo$receptor, ".", boo$sampletype),
  levels = paste0(rep(c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"), each = 2), ".", c("skin", "eAML"))
)


ggplot(boo[which(boo$sampletype %in% c("eAML", "skin")), ], aes(x = receptor.sampletype, y = value, color = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 0.5) +
  geom_signif(
    comparisons = list(
      c("TRA.skin", "TRA.eAML"),
      c("TRB.skin", "TRB.eAML"),
      c("TRD.skin", "TRD.eAML"),
      # c('TRG.skin', 'TRG.eAML'),
      c("IGH.skin", "IGH.eAML"),
      c("IGK.skin", "IGK.eAML"),
      c("IGL.skin", "IGL.eAML")
    ), color = "black", textsize = 3,
    test = "t.test",
    map_signif_level = function(x) {
      if (x < 0.001) {
        "<0.001"
      } else {
        as.character(round(x, digits = 3))
      }
    }
  ) +
  scale_y_sqrt("sequences", breaks = c(0, 10, 50, 100, 250, 500, 1000)) +
  scale_color_manual(values = c("eAML" = "orange", "skin" = "blue")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )
ggsave("./bulk_RNA_seq/figures/combined/plots/20241212_bulk_eAML_skin_immune_receptors.svg", width = 3, height = 2)
