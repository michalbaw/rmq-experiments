#!/bin/sh
make && python3 scripts/rmq_experiments.py --min_length=6 --max_length=7 --seq_type=random_walk --delta=0 --retries=100 #
Rscript ./scripts/my_rmq_stats.r #
Rscript ./scripts/cache_miss_stats.R #
