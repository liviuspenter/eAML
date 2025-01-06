library(dplyr)
library(ggplot2)
library(SummarizedExperiment)

Strobl.files <- data.table::fread("./data/Strobl/SraRunTable.txt")
source("./reanalysis_strobl/R/00_variant_calling_atac.R")

# create shell script for generating merged bams and calling mitochondrial DNA variants
fileConn <- file("./reanalysis_strobl/merge_bams.sh")
write("#!/bin/bash", file = "./reanalysis_strobl/merge_bams.sh", sep = "\n", append = T)

Strobl.metadata <- data.table::fread("./data/Strobl/SraRunTable.txt")
for (sample in unique(Strobl.metadata$BioSample)) {
  print(sample)
  bam.files <- paste0(
    Strobl.metadata$Run[which(Strobl.metadata$BioSample == sample)],
    ".mito/",
    Strobl.metadata$Run[which(Strobl.metadata$BioSample == sample)], ".mito.bam"
  )
  command <- paste0("samtools merge ", sample, ".merged.mito.bam ", paste(bam.files, collapse = " "))
  write(paste0("echo ", sample), file = "./reanalysis_strobl/merge_bams.sh", sep = "\n", append = T)
  write(command, file = "./reanalysis_strobl/merge_bams.sh", sep = "\n", append = T)
  write(paste0("mkdir ", sample), file = "./reanalysis_strobl/merge_bams.sh", sep = "\n", append = T)
  write(paste0("mv ", sample, ".merged.mito.bam ", sample), file = "./reanalysis_strobl/merge_bams.sh", sep = "\n", append = T)
  write(paste0("mgatk call -i ", sample, " -o ", sample, ".mgatk -c 8"), file = "./reanalysis_strobl/merge_bams.sh", sep = "\n", append = T)
}
close(fileConn)

# sample metadata
Strobl.metadata <- readxl::read_excel("./data/Strobl/samples.xlsx")
Strobl.metadata$label <- gsub(Strobl.metadata$label, pattern = "reanalysis blood_T_", replacement = "")
Strobl.metadata$label <- gsub(Strobl.metadata$label, pattern = "blood_T_", replacement = "")
Strobl.metadata$label <- gsub(Strobl.metadata$label, pattern = "skin_T_", replacement = "")
Strobl.metadata$patient <- stringr::str_split_fixed(Strobl.metadata$label, pattern = "_", n = 2)[, 1]
Strobl.metadata$day <- stringr::str_split_fixed(Strobl.metadata$label, pattern = "_", n = 2)[, 2]
Strobl.metadata$day <- gsub(Strobl.metadata$day, pattern = "day", replacement = "")
Strobl.metadata <- merge(Strobl.metadata, Strobl.files %>% group_by(`Sample Name`) %>% summarize(
  tissue.site = unique(tissue),
  bio.sample = unique(BioSample)
), by.x = "sample", by.y = "Sample Name")

# read variants
results.df <- data.frame()
for (s in unique(Strobl.metadata$bio.sample)) {
  message(s)
  results <- variant_calling_bulk(
    sample.bulk = readRDS(paste0("./data/Strobl/", s, ".mgatk/final/mgatk.rds")),
    coverage.position = 20, strand.coordination = 0, total.coverage = 30,
    sample.name = s
  )
  results <- results[which(!is.na(results$heteroplasmy)), ]
  results$sample <- s

  results.df <- rbind(results.df, results)
}

# manually inspect recipient-donor pairs with mismatched mtDNA variants
# 713, 718, 719, 721, 722, 728
# 712?

results.df <- merge(results.df, Strobl.metadata[, c("patient", "day", "tissue.site", "bio.sample")], by.x = "sample", by.y = "bio.sample")
results.df$tissue.site.day <- paste0(results.df$tissue.site, ".", results.df$day)
results.df$included <- ifelse(results.df$patient %in% c(713, 718, 719, 721, 722, 728), "yes", "no")

ggplot(results.df, aes(x = patient, y = coverage)) +
  geom_jitter(width = 0.2, size = 0.5) +
  geom_boxplot(aes(fill = included), outlier.size = 0, alpha = 0.5) +
  scale_y_continuous("mtDNA coverage") +
  scale_fill_manual(values = c("yes" = "blue", "no" = "grey90")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./reanalysis_strobl/figures/20241205_mtDNA_coverage_bulk.svg", width = 2.5, height = 2)

variants.all <- data.frame()

kinetics <- data.frame()
for (recipient in c(713, 718, 719, 721, 722, 728)) {
  message(recipient)
  boo <- as.data.frame(results.df[which(results.df$patient == recipient), c("variant", "tissue.site.day", "heteroplasmy")] %>%
    tidyr::pivot_wider(names_from = "tissue.site.day", values_from = "heteroplasmy", ))
  columns <- c("variant", "blood.-7", "blood.0", "skin.-7", "blood.14", "blood.100", "blood.365", "skin.14", "skin.100", "skin.365")
  columns <- columns[which(columns %in% colnames(boo))]
  boo <- boo[, columns]
  boo[is.na(boo)] <- 0
  rownames(boo) <- boo$variant
  donor.variants <- boo$variant[which(boo$blood.365 > 0.9 & boo$blood.0 < 0.1 & boo$`blood.-7` < 0.1 & boo$`skin.-7` < 0.1)]
  if (recipient != 722) {
    recipient.variants <- boo$variant[which(boo$blood.0 > 0.9 & boo$`blood.-7` > 0.9 & boo$`skin.-7` > 0.9 & boo$blood.365 < 0.1)]
  } else {
    recipient.variants <- boo$variant[which(boo$`skin.-7` > 0.9 & boo$`blood.-7` > 0.9 & boo$blood.365 < 0.1)]
  }

  boo$identity <- "none"
  boo[which(rownames(boo) %in% recipient.variants), "identity"] <- "recipient"
  boo[which(rownames(boo) %in% donor.variants), "identity"] <- "donor"

  boo <- boo[, -1]
  kinetics <- plyr::rbind.fill(kinetics, cbind(data.frame(recipient = recipient, variant = "donor"), t(data.frame(colMeans(boo[donor.variants, seq(1, 9)])))))
  kinetics <- plyr::rbind.fill(kinetics, cbind(data.frame(recipient = recipient, variant = "recipient"), t(data.frame(colMeans(boo[recipient.variants, seq(1, 9)])))))
  boo$patient <- recipient
  variants.all <- rbind(variants.all, boo)
}
# take out unplausible values
kinetics$blood.0[which(kinetics$recipient == 722 & kinetics$variant == "recipient")] <- NA
kinetics$blood.100[which(kinetics$recipient == 722 & kinetics$variant == "donor")] <- NA
kinetics.long <- kinetics %>% tidyr::pivot_longer(
  cols = c("blood.-7", "blood.0", "skin.-7", "blood.14", "blood.100", "blood.365", "skin.14", "skin.100", "skin.365"),
  names_to = "timepoint"
)
kinetics.long$recipient.variant <- paste0(kinetics.long$recipient, ".", kinetics.long$variant)
ggplot(
  kinetics.long[which(kinetics.long$timepoint %in% c("blood.-7", "blood.0", "blood.14", "blood.100", "blood.365")), ],
  aes(x = timepoint, y = 100 * value, color = variant)
) +
  geom_point(size = 0.5) +
  geom_line(aes(group = as.character(recipient.variant))) +
  scale_x_discrete("Day post-HSCT", limits = c("blood.-7", "blood.0", "blood.14", "blood.100", "blood.365"), labels = c("-7", "0", "14", "100", "365")) +
  scale_y_continuous("%mtDNA chimerism blood") +
  scale_color_manual(values = c("donor" = "purple", "recipient" = "orange")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./reanalysis_strobl/figures/20220824_Strobl_blood_longitudinal.svg", width = 2, height = 2)

ggplot(
  kinetics.long[which(kinetics.long$timepoint %in% c("skin.-7", "skin.14", "skin.100", "skin.365")), ],
  aes(x = timepoint, y = 100 * value, color = variant)
) +
  geom_point(size = 0.5) +
  geom_line(aes(group = as.character(recipient.variant))) +
  scale_x_discrete("Day post-HSCT", limits = c("skin.-7", NA, "skin.14", "skin.100", "skin.365"), labels = c("-7", "0", "14", "100", "365")) +
  scale_y_continuous("%mtDNA chimerism skin") +
  scale_color_manual(values = c("donor" = "purple", "recipient" = "orange")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./reanalysis_strobl/figures/20220824_Strobl_skin_longitudinal.svg", width = 2, height = 2)

# number of informative variants per case
stats <- variants.all %>%
  group_by(patient) %>%
  summarize(
    donor = length(which(identity == "donor")),
    recipient = length(which(identity == "recipient"))
  ) %>%
  tidyr::pivot_longer(cols = c("donor", "recipient"))

ggplot(stats, aes(x = as.character(patient), y = value, fill = name)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("donor" = "purple", "recipient" = "orange")) +
  scale_y_continuous("# mtDNA mutations") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./reanalysis_strobl/figures/20241205_number_of_variants_per_case.svg", width = 1.3, height = 2)

ggplot(variants.all[which(variants.all$identity != "none"), ], aes(x = identity, y = 100 * blood.14, color = identity)) +
  stat_summary(geom = "crossbar", width = 1, size = 0.5, fun.min = "median", fun.max = "median", fun = "median", color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_color_manual(values = c("donor" = "purple", "recipient" = "orange")) +
  scale_y_continuous("% heteroplasmy") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./reanalysis_strobl/figures/20241205_blood_chimerism_day14.svg", width = 1.1, height = 2)

ggplot(variants.all[which(variants.all$identity != "none"), ], aes(x = identity, y = 100 * skin.14, color = identity)) +
  stat_summary(geom = "crossbar", width = 1, size = 0.5, fun.min = "median", fun.max = "median", fun = "median", color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_color_manual(values = c("donor" = "purple", "recipient" = "orange")) +
  scale_y_continuous("% heteroplasmy") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./reanalysis_strobl/figures/20241205_skin_chimerism_day14.svg", width = 1.1, height = 2)

ggplot(variants.all[which(variants.all$identity != "none"), ], aes(x = identity, y = 100 * blood.365, color = identity)) +
  stat_summary(geom = "crossbar", width = 1, size = 0.5, fun.min = "median", fun.max = "median", fun = "median", color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_color_manual(values = c("donor" = "purple", "recipient" = "orange")) +
  scale_y_continuous("% heteroplasmy") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./reanalysis_strobl/figures/20241205_blood_chimerism_day365.svg", width = 1.1, height = 2)

ggplot(variants.all[which(variants.all$identity != "none"), ], aes(x = identity, y = 100 * skin.365, color = identity)) +
  stat_summary(geom = "crossbar", width = 1, size = 0.5, fun.min = "median", fun.max = "median", fun = "median", color = "black") +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_color_manual(values = c("donor" = "purple", "recipient" = "orange")) +
  scale_y_continuous("% heteroplasmy") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./reanalysis_strobl/figures/20241205_skin_chimerism_day365.svg", width = 1.1, height = 2)
