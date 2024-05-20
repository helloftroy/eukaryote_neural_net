#!/bin/bash
#SBATCH --job-name=star_genome
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_compare.csv
#SBATCH --error=job_output_compare.err

python compare_rice.py
