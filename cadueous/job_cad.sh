#!/bin/bash
#SBATCH --job-name=cad_ml
#SBATCH --partition=v100
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=c_out.txt
#SBATCH --error=c_err.txt

# run script
module load gcc/12.2.0
#gcc main.c -o main -std=c99
gcc *.c -o output -std=c99
gcc -std=c99
#~/.conda/envs/caduceus/bin/python -u train_caduceus.py
~/.conda/envs/caduceus/bin/python -u classify_caduceus.py

#nvidia-smi

