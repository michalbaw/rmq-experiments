# construction_times_plot.R

library(ggplot2)
library(grid)
library(plyr)
library(stringr) # <-- Add this line

# Custom theme required by the plotting function
theme_complete_bw <- function(base_size = 12, base_family = "") {
  theme(
    line =               element_line(colour = "black", size = 0.5, linetype = 1,
                                      lineend = "butt"),
    rect =               element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
    text =               element_text(family = base_family, face = "plain",
                                      colour = "black", size = base_size,
                                      hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                      margin = margin(), debug = FALSE),
    axis.text =          element_text(size = rel(0.8), colour = "grey50"),
    strip.text =         element_text(size = base_size * 0.7),
    axis.line =          element_blank(),
    axis.text.x =        element_text(size = base_size * 0.6 , lineheight = 0.9, angle = 0, colour = "black", vjust = 1),
    axis.text.y =        element_text(size = base_size * 0.7, lineheight = 0.9, colour = "black", hjust = 1),
    axis.ticks =         element_line(colour = "black"),
    axis.title.x =       element_blank(),
    axis.title.y =       element_text(size = base_size * 0.9, angle = 90, vjust = 0.5),
    axis.ticks.length =  unit(0.15, "cm"),
    
    legend.background =  element_blank(),
    legend.margin =      unit(0.25, "cm"),
    legend.key.height =  unit(0.5, "cm"),
    legend.key.width =   unit(0.5, "cm"),
    legend.text =        element_text(size = rel(0.75)),
    legend.text.align =  NULL,
    legend.title =       element_blank(),
    legend.title.align = NULL,
    legend.direction =   "horizontal",
    legend.justification = "center",
    legend.box =         NULL,
    legend.position =    "bottom",

    panel.background =   element_rect(fill = NA, colour = "grey", size = 1.3),
    panel.border =       element_blank(),
    panel.grid.major =   element_line(colour = "grey90", size = 0.7),
    panel.grid.minor =   element_line(colour = "grey90", size = 0.3),
    panel.margin =       unit(0.1, "lines"),

    strip.background =   element_rect(fill = NA, colour = NA),
    strip.text.x =       element_text(colour = "black", size = base_size * 0.8),
    strip.text.y =       element_text(colour = "black", size = base_size * 0.8, angle = -90),

    plot.background =    element_rect(colour = NA, fill = "white"),
    plot.title =         element_text(size = base_size * 1.2),
    plot.margin=         unit(c(3,3,3,3),"mm"),
    complete = TRUE
  )
}

# Plot which visualizes the construction time of the different algorithms
construction_time_plot <- function(c, title="") {

  c$ConstructTime <- as.numeric(as.character(c$ConstructTime))

  c$Algo <- revalue(c$Algo,
    c(
      "RMQ_FERRADA"="BP-Ferrada",
      "RMQ_SDSL_SCT"="SDSL-OLD",
      "RMQ_SUCCINCT"="SUCCINCT",
      "RMQ_SDSL_BP_FAST_REC_1024"="SDSL-BP-REC"
    )
  )
  algo_levels <- unique(c$Algo)

  algo_levels <- algo_levels[
    order(nchar(algo_levels), algo_levels)
  ]

  c$Algo <- factor(c$Algo, levels = algo_levels)

  # Aggregate statistics
  summary_df <- ddply(
    c,
    .(N, Algo),
    summarise,
    mean_time = mean(ConstructTime),
    sd_time = sd(ConstructTime),
    p95_time = quantile(ConstructTime, 0.95)
  )

  dodge <- position_dodge(width = 0.9)

  plot <- ggplot(
    summary_df,
    aes(
      x = factor(N),
      y = mean_time,
      fill = Algo
    )
  ) +
    ggtitle(title)

  # Mean bars
  plot <- plot +
    geom_bar(
      stat = "identity",
      position = dodge
    )

  # Standard deviation error bars
  plot <- plot +
    geom_errorbar(
      aes(
        ymin = mean_time - sd_time,
        ymax = mean_time + sd_time
      ),
      width = 0.2,
      position = dodge
    )

  # 95th percentile markers
  plot <- plot +
    geom_point(
      aes(
        y = p95_time,
        group = Algo
      ),
      position = dodge,
      size = 3,
      shape = 4
    )

  plot <- plot +
    scale_y_continuous(name = "Construction Time [ms]")

  plot <- plot +
    xlab("N")

  plot <- plot +
    theme_complete_bw()

  print(plot)
}
# ========== Data Loading and Execution =========== #

# Set your experiment directory path here
experiment_dir <- "./results/"
date <- "2026-06-08"
seq_type <- "random_walk"
max_length <- "6"
delta <- "0"
# tmp <- cbind(date,"rmq_experiment",seq_type,max_length,delta,"with_cache_misses")
tmp <- cbind(date,"rmq_experiment",seq_type,max_length,delta)
experiment <- str_c(tmp,collapse='_')
experiment_path <- paste(experiment_dir,experiment,sep="")

# Note: Update 'experiment_path' manually if the folder structure is different
# experiment_path <- "path/to/your/experiment/folder"

# Read the construction results data
c <- read.csv2(paste(experiment_path, "/construct_result.csv", sep=""), sep=",", header=TRUE)

# Filter out specific algorithms as done in the original script
c <- subset(c, c$Algo != "RMQ_SDSL_BP_FAST_REC_1024")

# Format numeric columns
c$N <- as.numeric(as.character(c$N))
c$BPE <- as.numeric(as.character(c$BPE))

# Generate the plot (Uncommented from your original script)
construction_time_plot(c, title="Algorithm Construction Times")