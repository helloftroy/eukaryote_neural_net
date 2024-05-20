#!/bin/bash
#SBATCH --job-name=aminoacid
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err

# usage: the (1) intron.fasta file and (2) the .igv.con.intron.bed file
awk -F',' '{print $1 "," $6}' ${1} > extracted_columns.csv
awk '{print $1":"$2"-"$3"("$6")", $4}' ${2} > new_formatted_file.txt

sort -t',' -k1,1 extracted_columns.csv > sorted_extracted_columns.csv
awk -F' ' -v OFS=',' '{print $2, $1}' new_formatted_file.txt | sort > sorted_new_formatted_file.txt
join -t',' -1 1 -2 2 -o 1.1,1.2,2.1,2.2 sorted_extracted_columns.csv sorted_new_formatted_file.txt > merged_file.csv
sort -t',' -k1,1 sorted_new_formatted_file.txt > sorted_new_formatted_file_sorted.txt  
join -t',' -1 1 -2 1 -o 1.1,1.2,2.2 sorted_extracted_columns.csv sorted_new_formatted_file_sorted.txt > merged_file.csv

awk -F',' 'BEGIN {print "intron_id,sa.eff"} {print $3 "," $2}' merged_file.csv > maturation_eff.csv



#awk -F',' '{print $1 "," $7}' intron.info.csv > genes_list.txt
#sort -t ',' -k1,1 -V genes_list.txt > sorted_genes_list.txt

#awk 'BEGIN{FS=","} NR==FNR{genes[$1]; next} $1 ~ /^>/ {match($1, />(.*)/, m); if (m[1] in genes) {p=1} else {p=0}} p' sorted_genes_list.txt modified_headers.fasta > genes_of_interest.fasta

