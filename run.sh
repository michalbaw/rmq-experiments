#!/bin/bash
make && python3 scripts/rmq_experiments.py --min_length=5 --max_length=7 --seq_type=random --delta=0 --count_cache_misses=1 --retries=10
Rscript ./scripts/my_rmq_stats.r
Rscript ./scripts/cache_miss_stats.R
