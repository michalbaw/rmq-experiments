library(dplyr)
library(ggplot2)

data <- read.csv("results/rmq_experiment/query_result.csv")

stats <- data %>%
  group_by(Algo, N, Range) %>%
  summarise(
    count = n(),
    mean_time = mean(Time),
    sd_time = sd(Time),
    min_time = min(Time),
    max_time = max(Time),
    q25 = quantile(Time, 0.25),
    median = quantile(Time, 0.50),
    q75 = quantile(Time, 0.75),
    q95 = quantile(Time, 0.95),
    q99 = quantile(Time, 0.99),
    .groups = 'drop'
  )

# print(as.data.frame(stats), row.names = FALSE)

write.csv(stats, "rmq_statistics.csv", row.names = FALSE)

cat("\n=== Summary by Algorithm ===\n")
for (algo in unique(data$Algo)) {
  cat("\n", algo, ":\n", sep = "")
  algo_stats <- stats %>% filter(Algo == algo)
  # print(as.data.frame(algo_stats), row.names = FALSE)
}

overall_stats <- data %>%
  group_by(Algo) %>%
  summarise(
    count = n(),
    mean_time = mean(Time),
    sd_time = sd(Time),
    median = median(Time),
    .groups = 'drop'
  )

cat("\n=== Overall Statistics by Algorithm ===\n")
# print(as.data.frame(overall_stats), row.names = FALSE)

cat("\n=== Average Times for N = 10^6, Range = 10^5 ===\n")
specific_stats <- stats %>%
  filter(N %in% c(1000000, 1e6), Range %in% c(100000, 1e5)) %>%
  select(Algo, mean_time)
print(as.data.frame(specific_stats), row.names = FALSE)

cat("\n=== Generating Plots ===\n")
n_values <- unique(stats$N)
plot_list <- list()

for (n_val in n_values) {
  cat("Creating plot for N =", n_val, "\n")
  plot_data <- stats %>% filter(N == n_val)

  plot_data <- stats %>%
    filter(N == n_val) %>%
    filter(Algo %in% c(
      # "RMQ_SDSL_SCT",
      # "RMQ_SUCCINCT",
      # "RMQ_FERRADA",
      # "RMQ_SDSL_REC"
      "RMQ_SDSL_FAST"
      ,"RMQ_SDSL_REC_ST"
      # ,"RMQ_SDSL_FAST_ST"
      ,"RMQ_SDSL_SPARSE_BITMASKS"
    ))
  
  p <- ggplot(plot_data, aes(x = Range, y = mean_time, color = Algo, group = Algo)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_log10(labels = scales::scientific) +
    scale_y_continuous(limits = c(0, 0.35), expand = expansion(mult = c(0, 0.05))) +
    # scale_y_continuous() + 
    scale_color_manual(
      values = c(
        "RMQ_SDSL_SCT" = "#E69F00",
        "RMQ_SUCCINCT" = "#56B4E9",
        "RMQ_FERRADA" = "#009E73",
        "RMQ_SDSL_REC" = "#F0E442",
        "RMQ_SDSL_FAST" = "#0072B2",
        "RMQ_SDSL_REC_ST" = "#D55E00",
        "RMQ_SDSL_FAST_ST" = "#CC79A7",
        "RMQ_SDSL_SPARSE_BITMASKS" = "#9600d5"
      ),
      breaks = c(
        "RMQ_SDSL_SCT",
        "RMQ_SUCCINCT",
        "RMQ_FERRADA",
        "RMQ_SDSL_REC",
        "RMQ_SDSL_FAST",
        "RMQ_SDSL_REC_ST",
        "RMQ_SDSL_FAST_ST",
        "RMQ_SDSL_SPARSE_BITMASKS"
      ),
      labels = c(
        "RMQ_SDSL_FAST" = "Alstrup et al.",
        "RMQ_SDSL_FAST_ST" = "Alstrup et al. + sparse table",
        "RMQ_SDSL_REC" = "Baumstark et al.",
        "RMQ_SDSL_REC_ST" = "Baumstark et al.",
        "RMQ_SUCCINCT" = "Succinct (memory-optimal)",
        "RMQ_FERRADA" = "Ferrada & Navarro",
        "RMQ_SDSL_SCT" = "SDSL library",
        "RMQ_SDSL_SPARSE_BITMASKS" = "Hybrid Baumstark + Alstrup"
      )
    ) +
    labs(
      title = bquote("RMQ Algorithm Performance (N =" ~ .(scales::scientific(n_val)) ~ ")"),
      x = "Query Range",
      y = "Mean Time (microseconds)",
      color = "Algorithm"
    ) +
    theme_minimal(base_size = 16) +
    # theme(
    #   plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    #   axis.title = element_text(size = 18),
    #   axis.text = element_text(size = 14),
    #   legend.title = element_text(size = 16),
    #   legend.text = element_text(size = 14),
    #   legend.position = "right"
    # )
    theme(
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      axis.title.x = element_text(size = 18,  margin = margin(t = 10)),
      axis.title.y = element_text(size = 18, margin = margin(r = 10)),
      
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.position = "right",
      
      legend.key.size = unit(1, "cm")
    )

  plot_list[[length(plot_list) + 1]] <- p
  
  # Save individual plot
#   filename <- paste0("rmq_plot_N_", n_val, ".png")
#   ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)
}

# Save all plots in a single PDF file
cat("\nSaving all plots to single PDF file...\n")
pdf("rmq_all_plots.pdf", width = 10, height = 6)
for (p in plot_list) {
  print(p)
}
dev.off()
cat("All plots saved together in: rmq_all_plots.pdf\n")

construct_data <- read.csv("results/rmq_experiment/construct_result.csv")

cat("=== Construction Data Summary ===\n")
print(summary(construct_data))

construct_stats <- construct_data %>%
  group_by(Algo, N) %>%
  summarise(
    count = n(),
    mean_construct_time = mean(ConstructTime),
    sd_construct_time = sd(ConstructTime),
    mean_bpe = mean(BPE),
    sd_bpe = sd(BPE),
    .groups = 'drop'
  )

cat("\n=== Construction Statistics ===\n")
print(as.data.frame(construct_stats), row.names = FALSE)

cat("\n=== Creating BPE Plot ===\n")

bpe_plot <- ggplot(construct_stats, aes(x = N, y = mean_bpe, color = Algo, group = Algo)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_x_log10(labels = scales::scientific) +
  scale_y_continuous() +
  labs(
    title = "RMQ Algorithm Memory Usage (Bits Per Element)",
    x = "N (Structure Size)",
    y = "Bits Per Element (BPE)",
    color = "Algorithm"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  )

# Save BPE plot
ggsave("rmq_bpe_plot.pdf", plot = bpe_plot, width = 10, height = 6)
