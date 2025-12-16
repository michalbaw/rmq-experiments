#!/bin/bash
make && python3 scripts/rmq_experiments.py --min_length=5 --max_length=8 --seq_type=random --delta=0