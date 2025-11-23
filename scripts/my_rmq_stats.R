# Install and load required libraries
install.packages("ggplot2", repos = "http://cran.us.r-project.org")
library(dplyr)
library(ggplot2)

# Read the CSV file
data <- read.csv("results/2025-11-23_rmq_experiments_random_8_0.csv")

# Calculate statistics for each algorithm and N
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

# Print the results (force full display)
print(as.data.frame(stats), row.names = FALSE)

# Optional: Save results to CSV
write.csv(stats, "rmq_statistics.csv", row.names = FALSE)

# Optional: Create a summary table for each algorithm
cat("\n=== Summary by Algorithm ===\n")
for (algo in unique(data$Algo)) {
  cat("\n", algo, ":\n", sep = "")
  algo_stats <- stats %>% filter(Algo == algo)
  print(as.data.frame(algo_stats), row.names = FALSE)
}

# Optional: Overall statistics by algorithm (across all N)
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

# Generate plots for each N value
cat("\n=== Generating Plots ===\n")

# Get unique N values
n_values <- unique(stats$N)

# Create a list to store all plots
plot_list <- list()

# Create a plot for each N value
for (n_val in n_values) {
  cat("Creating plot for N =", n_val, "\n")
  
  # Filter data for this N
  plot_data <- stats %>% filter(N == n_val)
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Range, y = mean_time, color = Algo, group = Algo)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_continuous(labels = scales::scientific) +
    scale_y_continuous(labels = scales::scientific) +
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
  
  # Add to list
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