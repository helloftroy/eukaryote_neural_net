#!/bin/bash
#SBATCH --job-name=final_score
#SBATCH --cpus-per-task=1
#SBATCH --mem=9G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err

#usage qsub S8_intron.score.sh $genome.dir $intron_info_csv $maturation.eff.csv
# for the intron info, only need to use the TT info

genome="${1%/}"

intronbed="${genome}/*all.intron.bed"
intronbytxbed="${genome}/*intronbytx.bed"
e10i12bed="${genome}/*e10i12.bed"
e10i12seq="${genome}/*e10i12.seq"

module load R/4.2.2
Rscript  ~/tools/S8_intron.score.R $2 $3 $intronbed $intronbytxbed $e10i12bed $e10i12seq

