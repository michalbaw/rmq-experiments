SDSL RMQ Experiment Repository
=======

Install
-----
```
git clone https://github.com/kittobi1992/rmq-experiments.git
git submodule init
git submodule update
make
```

Execute RMQ-Experiments
---
Make sure you have created the folders results/ and benchmark/ before you start the experiment script. 

Starting the experiment:
```
python3 scripts/rmq_experiments.py --min_length=5 --max_length=8 --seq_type=random --delta=0
```
* ``--length=N`` specifies the size of the sequence, which is 10<sup>N</sup>
* ``--seq_type=type`` with ``type = {random,increasing,decreasing}`` specifies the input type (see paper)
* ``--delta=x`` specifies the delta of pseudo-increasing and -decreasing sequences (see paper) 

Results are located in folder results/ after experiment execution. 
Use R-Script ``script/rmq_stats.R`` to convert results into tikz-Figures.

Execute LCP-Experiments
---
Make sure you have downloaded the pizza&chilli benchmarks with script ``pizza&chilli/download.sh``

Starting the experiment:
```
python scripts/lcp_experiments.py
```

Results are located in folder results/ after experiment execution. 
Use R-Script ``script/lcp_stats.R`` to convert results into tikz-Figures.
