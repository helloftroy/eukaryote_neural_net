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

python extract_tss_tts.py

module load bedops/2.4.39
module load bedtools/2-2.29.2

#awk '$4 >= 0' tss_tts_output.gff > filtered_tss_tts_output.gff

#gff2bed < filtered_tss_tts_output.gff > utr_output.bed
bedtools getfasta -fi ${1} -bed utr_output.bed -s -fullHeader -fo tss_sequences.fa 
bedtools getfasta -fi ${1} -bed utr_output.bed -s -tab -fo tss_sequences.csv

paste <(awk '{print $9, $3}' tss_tts_output.gff) tts.csv | tr ' ' '\t' > tts.csv

paste -d '\t' tts.csv tss_sequences.csv > tss_tts_sequences_with_names.csv
