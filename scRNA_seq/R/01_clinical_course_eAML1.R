WBC.data <- as.data.frame(readxl::read_excel("./data/clinical_data/20230118_eAML1.xlsx"))

p <- ggplot() +
  geom_line(data = WBC.data, aes(x = Day, y = WBC)) +
  geom_point(data = WBC.data, aes(x = Day, y = WBC), size = 0.5) +
  scale_x_continuous("Day after initial diagnosis") +
  scale_y_continuous("WBC [1,000 / µl]", limits = c(0, 60)) +
  geom_point(data = WBC.data[which(!is.na(WBC.data$Label)), ], aes(x = Day), y = 40, size = 1, shape = 25) +
  annotate("text",
    x = WBC.data[which(!is.na(WBC.data$Label)), "Day"], y = 50,
    label = WBC.data[which(!is.na(WBC.data$Label)), "Label"], angle = 90, size = 2
  ) +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./scRNA_seq/figures/clinical_plots/20230118_eAML1_WBC.svg", width = 2.5, height = 2, plot = p)
