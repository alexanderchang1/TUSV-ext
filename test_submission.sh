python tusv-ext.py \
  -i /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext/simulation_data/experiment_3_1_100_20_n/patient1/sample \
  -o /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext/simulation_data/experiment_3_1_100_20_n/patient1/output \
  -n 5 \
  -c 10 \
  -t 1000 \
  -r 10 \
  -col \
  -leaf 



python tusv-ext.py \
  -i /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext/simulation_data/experiment_3_1_100_20_n/patient1/sample \
  -o /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext/simulation_data/experiment_3_1_100_20_n/patient1/sample/output \
  -n 10 \  # Set to a reasonable maximum number of nodes
  -c 10 \  # Example value for maximum copy number
  -t 1000 \  # Example value for maximum coordinate-descent iterations
  -r 10 \  # Example value for random initializations
  -col \  # Enable collapsing redundant nodes
  -b  # Automatically set hyper-parameters lambda 1 and lambda 2 (recommended)

# Single cell sequencing of metastatic lesions have shown various leiden clusters, so recommend 5 x tissue sample, so for example if you have 4 tissue samples n = 20.