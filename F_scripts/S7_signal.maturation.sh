#!/bin/bash
#SBATCH --job-name=signal_mature
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err

#usage qsub S7_signal.maturation.sh $genome.dir output.name
#this will auto detect all txbyexon.signal.dt files in current folder and divide them into rna or tt. copy paired signal files to current folder!

genome="${1%/}"

exonbytxfile="${genome}/*exonbytx.bed"

module load R/4.2.2
Rscript  ~/tools/S7_signal.maturation.R $exonbytxfile $2

