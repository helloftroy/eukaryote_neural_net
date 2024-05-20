#!/bin/bash
#SBATCH --job-name=star_genome
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err



STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir . \
--genomeFastaFiles $1 \
--sjdbGTFfile $2 \
--sjdbOverhang 149 \
--sjdbGTFtagExonParentTranscript Parent
