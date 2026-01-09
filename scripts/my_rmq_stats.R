library(dplyr)
library(ggplot2)

data <- read.csv("results/2026-01-09_rmq_experiment_random_6_0_with_cache_misses/query_result.csv")

data <- data %>%
  mutate(RangeBin = 10^floor(log10(Range)))

stats <- data %>%
  group_by(Algo, N, RangeBin) %>%
  summarise(
    count = n(),
    mean_time = mean(Time),
    sd_time = sd(Time),
    sem_time = sd_time / sqrt(count),
    ci_low = mean_time - 1.96 * sem_time,
    ci_high = mean_time + 1.96 * sem_time,
    min_time = min(Time),
    max_time = max(Time),
    q025 = quantile(Time, 0.025),
    q975 = quantile(Time, 0.975),
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
    q025 = quantile(Time, 0.025),
    q975 = quantile(Time, 0.975),
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
  
  plot_data <- stats %>%
    filter(N == n_val) %>%
    filter(!Algo %in% c("RMQ_SDSL_SCT", "RMQ_SUCCINT", "RMQ_FAST", "RMQ_FERRADA"))

  p <- ggplot(plot_data, aes(x = RangeBin, y = mean_time, color = Algo, group = Algo)) +
  geom_ribbon(
    aes(
      ymin = q025,
      ymax = q975,
      fill = Algo
    ),
    alpha = 0.2,
    color = NA
  ) + 
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.1,
    alpha = 0.6
  ) +
  scale_x_log10(labels = scales::scientific) +
  labs(
    title = bquote("RMQ Algorithm Performance (N =" ~ .(scales::scientific(n_val)) ~ ")"),
    x = "Range Bin [10^k, 10^(k+1))",
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
