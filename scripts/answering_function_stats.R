path_base <- "results/2026-01-13_rmq_experiment_random_8_0/"

data <- read.csv(paste(path_base, "query_result.csv", sep=""))

filtered_data_rec <- data[data$Algo == "RMQ_SDSL_REC", ]
filtered_data_rec_st <- data[data$Algo == "RMQ_SDSL_REC_ST", ]

library(tidyr)
library(dplyr)

get_mode_and_pct <- function(x) {
  if (length(x) == 0) return(NA)
  
  counts <- table(x)
  mode_val <- names(counts)[which.max(counts)]
  pct <- max(counts) / length(x) * 100

  return(sprintf("%s (%.1f%%)", mode_val, pct))
}

# RMQ_SDSL_REC
result_table_rec <- filtered_data_rec %>%
  group_by(Range, N) %>%
  summarise(Function_Mode = get_mode_and_pct(Function), .groups = "drop") %>%
  pivot_wider(
    id_cols = Range,
    names_from = N,
    values_from = Function_Mode
  )

# RMQ_SDSL_REC_ST
result_table_rec_st <- filtered_data_rec_st %>%
  group_by(Range, N) %>%
  summarise(Function_Mode = get_mode_and_pct(Function), .groups = "drop") %>%
  pivot_wider(
    id_cols = Range,
    names_from = N,
    values_from = Function_Mode
  )

cat("\n=== RMQ_SDSL_REC ===\n")
print(result_table_rec, width = Inf)

cat("\n\n=== RMQ_SDSL_REC_ST ===\n")
print(result_table_rec_st, width = Inf)