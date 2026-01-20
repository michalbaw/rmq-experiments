#!/bin/bash
make && python3 scripts/rmq_experiments.py --min_length=2 --max_length=8 --seq_type=random --delta=0 --count_cache_misses=0 --check_function=1