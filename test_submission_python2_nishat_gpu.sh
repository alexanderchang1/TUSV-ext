#!/bin/bash

#SBATCH --job-name=TUSV_EXT_PYTHON2
#SBATCH --output=TUSV_EXT_PYTHON2_%j.log
#SBATCH --error=TUSV_EXT_PYTHON2_%j.log
#SBATCH --time=12:00:00
#SBATCH --cluster=gpu
#SBATCH --partition=a100,a100_nvlink
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL         # Email notifications for job done/failed
#SBATCH --mail-user=alc376@pitt.edu  # Replace with your email



echo "Starting job $SLURM_JOB_ID at $(date)"

source /bgfs/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/tusv_ext_python2.7/
echo "Activated conda environment"

module load gurobi/11.0.2

echo "Loaded gurobi module"

# Single cell sequencing of metastatic lesions have shown various leiden clusters, so recommend 5 x tissue sample, so for example if you have 4 tissue samples n = 20.

cd /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext
echo "Changed directory to /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext"

python tusv-ext.py \
  -i /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext/aer5_frozen/AER5/ \
  -o /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/TUSV-ext/aer5_frozen/AER5_output/ \
  -n 3 -c 10 -t 2 -r 2 -m 1000 -col -b -sv_ub 80 -C 120
echo "Finished running tusv-ext.py"

echo "Job $SLURM_JOB_ID completed at $(date)"