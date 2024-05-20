#!/bin/bash
#SBATCH --job-name=node_test
#SBATCH --partition=v100
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9G
#SBATCH --output=out.txt
#SBATCH --error=er.txt

# run script
#~/.conda/envs/caduceus_env/bin/python -u train_caduceus.py
~/.conda/envs/corn/bin/python -u classifier_withstart.py
#~/.conda/envs/corn/bin/python -u train.py
