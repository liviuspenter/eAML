# generate swimers plots and oncoplot

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

source("./scRNA_seq/00_eAML.colors.R")

metadata <- readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 4)
metadata$diagnosis.timepoint <- as.numeric(metadata$eAML_relapse - metadata$Diagnosis)
metadata$transplant.timepoint <- as.numeric(metadata$eAML_relapse - metadata$Transplant)
metadata$death.timepoint <- as.numeric(metadata$Death - metadata$eAML_relapse)
metadata$death.timepoint[is.na(metadata$Death)] <- as.numeric(as.POSIXct("2024-08-01") - metadata$eAML_relapse[is.na(metadata$Death)]) # use 08/01/2024 for patients still alive

samples <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 5))
for (i in seq(1, nrow(samples))) {
  samples[i, "sample.timepoint"] <- as.numeric(as.POSIXct(samples[i, "Date"]) - as.POSIXct(metadata$eAML_relapse[which(metadata$eAML == samples$eAML[i])]))
}

events <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 6))
for (i in seq(1, nrow(events))) {
  events[i, "event.timepoint"] <- difftime(as.POSIXct(events[i, "Date"]), as.POSIXct(metadata$eAML_relapse[which(metadata$eAML == events$eAML[i])]))
}

eAML.order <- order(metadata$death.timepoint)

metadata <- metadata[eAML.order, ]
metadata$rank <- as.numeric(seq(1, nrow(metadata)))
samples$rank <- as.numeric(sapply(samples$eAML, FUN = function(x) {
  metadata[which(metadata$eAML == x), "rank"]
}))
events$rank <- as.numeric(sapply(events$eAML, FUN = function(x) {
  metadata[which(metadata$eAML == x), "rank"]
}))

# swimmers plot

ggplot() +
  geom_segment(data = metadata, aes(x = rank, xend = rank, y = -diagnosis.timepoint, yend = 0), color = "darkgrey") +
  geom_segment(data = metadata[which(!is.na(metadata$Transplant)), ], aes(x = rank, xend = rank, y = -transplant.timepoint, yend = 0), color = "#00aeef") +
  geom_segment(data = metadata, aes(x = rank, xend = rank, y = 0, yend = death.timepoint), color = "firebrick") +
  geom_point(data = metadata[which(!is.na(metadata$Death)), ], aes(x = rank, y = death.timepoint + 100), color = "black", shape = 3) +
  geom_point(data = samples, aes(x = rank - 0.1, y = sample.timepoint, fill = tissue), shape = 21, color = "black") +
  geom_point(data = events, aes(x = rank + 0.1, y = event.timepoint, shape = event)) +
  scale_x_continuous(breaks = seq(1, 7), labels = paste0("eAML", eAML.order)) +
  scale_y_continuous("Days from eAML diagnosis") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  scale_shape_manual(values = c("chemotherapy" = 25, "DLI" = 17, "DLI_ipilimumab" = 15, "ipilimumab" = 18, "nivolumab" = 20, "radiotherapy" = 8, "transplant" = 9)) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.y = element_blank()
  )
ggsave("./scRNA_seq/figures/clinical_plots/20240802_eAML_swimmers_plot.svg", width = 7, height = 2)

ggplot() +
  geom_segment(data = metadata, aes(x = rank, xend = rank, y = -diagnosis.timepoint, yend = 0), color = "darkgrey") +
  geom_segment(data = metadata[which(!is.na(metadata$Transplant)), ], aes(x = rank, xend = rank, y = -transplant.timepoint, yend = 0), color = "#00aeef") +
  geom_segment(data = metadata, aes(x = rank, xend = rank, y = 0, yend = death.timepoint), color = "firebrick") +
  geom_point(data = metadata[which(!is.na(metadata$Death)), ], aes(x = rank, y = death.timepoint + 100), color = "black", shape = 3) +
  geom_point(data = samples, aes(x = rank - 0.1, y = sample.timepoint, fill = tissue), shape = 21, color = "black") +
  geom_point(data = events, aes(x = rank, y = event.timepoint), shape = 124) +
  scale_x_continuous(breaks = seq(1, 7), labels = paste0("eAML", eAML.order)) +
  scale_y_continuous("days from eAML") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  scale_shape_manual(values = c("chemotherapy" = 25, "DLI" = 17, "DLI_ipilimumab" = 15, "ipilimumab" = 18, "nivolumab" = 20, "radiotherapy" = 8, "transplant" = 9)) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.y = element_blank()
  )
ggsave("./scRNA_seq/figures/clinical_plots/20241230_eAML_swimmers_plot.svg", width = 4, height = 2)

### oncoplot

oncoplot.data <- as.data.frame(readxl::read_excel("./data/20240802_metadata_eAML.xlsx", sheet = 7))
oncoplot.data$dummy <- 1
oncoplot.data <- as.data.frame(oncoplot.data %>% tidyr::pivot_wider(names_from = "Gene", values_from = "dummy"))
rownames(oncoplot.data) <- oncoplot.data$eAML
oncoplot.data <- oncoplot.data[, -c(1, 2, 3)]
oncoplot.data <- oncoplot.data[, order(colSums(oncoplot.data, na.rm = T), decreasing = T)]
oncoplot.data[is.na(oncoplot.data)] <- 0
oncoplot.data <- t(oncoplot.data)

col_fun <- circlize::colorRamp2(breaks = c(0, 1), colors = c("grey90", "firebrick"))
svglite::svglite("./scRNA_seq/figures/clinical_plots/20240802_eAML_oncoplot.svg", width = 1.6, height = 2.3)
Heatmap(oncoplot.data,
  cluster_columns = F, cluster_rows = F, row_names_side = "left",
  column_labels = paste0("eAML", colnames(oncoplot.data)), border = T, col = col_fun, show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10)
)
dev.off()
