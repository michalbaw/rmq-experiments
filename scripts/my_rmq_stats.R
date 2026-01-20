library(dplyr)
library(ggplot2)

path_base <- "results/2025-12-17_rmq_experiment_random_8_0/"

data <- read.csv(paste(path_base, "query_result.csv", sep=""))

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

print(as.data.frame(stats), row.names = FALSE)

write.csv(stats, "rmq_statistics.csv", row.names = FALSE)

cat("\n=== Summary by Algorithm ===\n")
for (algo in unique(data$Algo)) {
  cat("\n", algo, ":\n", sep = "")
  algo_stats <- stats %>% filter(Algo == algo)
  print(as.data.frame(algo_stats), row.names = FALSE)
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
print(as.data.frame(overall_stats), row.names = FALSE)

cat("\n=== Generating Plots ===\n")
n_values <- unique(stats$N)
plot_list <- list()

for (n_val in n_values) {
  cat("Creating plot for N =", n_val, "\n")
  plot_data <- stats %>% filter(N == n_val)

  plot_data <- stats %>%
    filter(N == n_val) %>%
    filter(!Algo %in% c("RMQ_SDSL_SCT", "RMQ_SUCCINT", "RMQ_FAST", "RMQ_FERRADA"))
  
  p <- ggplot(plot_data, aes(x = Range, y = mean_time, color = Algo, group = Algo)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_log10(labels = scales::scientific) +
    scale_y_continuous() +
    labs(
      title = bquote("RMQ Algorithm Performance (N =" ~ .(scales::scientific(n_val)) ~ ")"),
      x = "Range (Query Range)",
      y = "Mean Time (seconds)",
      color = "Algorithm"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
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

cat("\nIndividual plots saved as: rmq_plot_N_*.png\n")
cat("All plots saved together in: rmq_all_plots.pdf\n")

construct_data <- read.csv(paste(path_base, "construct_result.csv", sep=""))

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