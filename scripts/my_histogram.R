library(dplyr)
library(ggplot2)

data <- read.csv("results/2026-01-09_rmq_experiment_random_7_0_with_cache_misses/query_result.csv")

filtered_data <- data %>%
  filter(
    N == 1e5,
    Algo == "RMQ_SDSL_FAST",
    Range == 1e4
  )

stopifnot(nrow(filtered_data) > 0)

p95 <- quantile(filtered_data$Time, 0.95)

filtered_p95 <- filtered_data %>%
  filter(Time <= p95)

p <- ggplot(filtered_p95, aes(Time)) +
  geom_histogram(bins = 25) +
  coord_cartesian(xlim = c(0, 0.25)) +
  labs(
    title = "Histogram of Query Time (≤ 95th Percentile)",
    subtitle = sprintf(
      "N = 1e5, Algo = RMQ_SDSL_FAST, Range = 1e4\n95th percentile = %.3g",
      p95
    ),
    x = "Time (seconds)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  filename = "rmq_histogram_p95_x_0_1_N_1e5_RMQ_SDSL_FAST_Range_1e4.pdf",
  plot = p,
  width = 8,
  height = 6
)

cat("PDF saved successfully using ggsave().\n")
