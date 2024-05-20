#!/bin/bash
#SBATCH --job-name=aminoacid
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output.out
#SBATCH --error=job_output_.err


#grep -E 'gene|polypeptide' rice_japonica_v3.igv.gff > filtered.gff
#sed -E 's/(ID=[^;-]+)-[0-9]+;/\1;/g' filtered.gff > filtered_output.gff

python faster_extract_3UTR_5UTR.py

module load bedops/2.4.39
module load bedtools/2-2.29.2

gff2bed < 53utr_output.gff > 53utr_output.bed
bedtools getfasta -fi ${1} -bed 53utr_output.bed -s -fullHeader -fo utr_sequences.fa 
bedtools getfasta -fi ${1} -bed 53utr_output.bed -s -tab -fo utr_sequences.csv

paste <(awk '{print $9, $3}' 53utr_output.gff) genes.csv | tr ' ' '\t' > genes.csv
awk -F';' '{print $2,$3}' genes.csv > new_genes.csv

paste -d '\t' new_genes.csv utr_sequences.csv > utr_sequences_with_names.csv
