#!/bin/bash
#SBATCH --job-name=aminoacid
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_j.out
#SBATCH --error=job_output_j.err

# just take out the first intron from the intron file
#grep -E -A 1 "intron1($|[^0-9])" output.fasta > intron_one.fasta

# Extractithe specified columns from the original file and save them to a new file
cut -d',' -f2,7,9 rice.maturation_eff_exp_genes.csv > extracted_columns.csv
cut -d',' -f2,7,9 batch5_rice.maturation_eff_exp_genes.csv > 5extracted_columns.csv

# Sort the columns based on the values in descending order
sort -t',' -k2,2nr -k3,3nr extracted_columns.csv > sorted_extracted_columns.csv

# Calculate the number of lines in the file and top 10%
total_lines=$(wc -l < sorted_extracted_columns.csv)
top_percentage_lines=$((total_lines * 6 / 100))

# Select the top 10% of the sorted data into a new file
head -n "$top_percentage_lines" sorted_extracted_columns.csv > top_10_percent_sorted_extracted_column.csv
awk -F',' '$3 >= 1.07' top_10_percent_sorted_extracted_column.csv > filtered_top_10_percent_sorted_extracted_column.csv
cut -d ',' -f 1 filtered_top_10_percent_sorted_extracted_column.csv > genes_list.txt

# Filter the rows of 5extracted_columns.csv based on genes_list.txt
awk -F ',' 'FNR==NR{genes[$1]; next} $1 in genes' genes_list.txt 5extracted_columns.csv > filtered_5.csv
awk -F ',' '($3 < 0 || $2 < 0.05) {print $1}' filtered_5.csv > ids_to_remove.txt

# filter out the genes that are < 0 maturation factor in the other batch
awk -F ',' 'FNR==NR{ids[$1]; next} !($1 in ids)' ids_to_remove.txt filtered_top_10_percent_sorted_extracted_column.csv > updated.csv
#pull out the intron 1 from the fasta file
#sed 's/$/.intron1/' genes_list.txt > mod_genes_list.txt
#awk 'BEGIN{while((getline line < "mod_genes_list.txt") > 0) patterns[line]=1} /^>/{p=0} {if(substr($1,2) in patterns)p=1} p' intron_one.fasta > intron_one_me.fasta
sed -i 's/t/g/g' updated.csv
awk -F ',' '{print $1}' updated.csv | grep -w -F -f - utr_sequences_with_names.csv > utr.csv
awk 'BEGIN {FS="\t"; OFS="\t"} {len=length($5); printf "%s\t%s\t%d\t%s\n", $1, $2, len, $5}' utr.csv > length_top_utr.csv
awk -F'\t' '$3 >= 19' length_top_utr.csv > filtered_length_top_utr.csv

wc -l filtered_length_top_utr.csv
awk -F'\t' '{sum += $3} END {print sum}' filtered_length_top_utr.csv

rm genes_list.txt 5extracted_columns.csv filtered_5.csv 5extracted_columns.csv length_top_utr.csv filtered_top_10_percent_sorted_extracted_column.csv top_10_percent_sorted_extracted_column.csv mod_genes_list.txt mod_two_genes_list.txt sorted_extracted_columns.csv two_sorted_extracted_columns.csv two_genes_list.txt two_sorted_extracted_columns.csv sorted_extracted_columns.csv extracted_columns.csv

 
