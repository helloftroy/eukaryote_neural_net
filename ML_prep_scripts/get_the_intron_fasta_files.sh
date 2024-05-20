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
#awk '{print $1":"$2"-"$3"("$6")", $4}' ${2} > new_formatted_file.txt
awk 'NR==FNR{a[">"$1]=$2;next}{if($0~/>/){sub(/>.*/,">"a[substr($0,2)]);print}else{printf "%s",$0}}' new_formatted_file.txt ${1} > modified_headers.fasta
#awk 'NR==FNR{a[">"$1]=$2;next}{if($0~/>/){sub(/>.*/,">"a[$1]);print}else{print}}' new_formatted_file.txt ${1} > modified_headers.fasta
#awk -F',' '{print $1 "," $7}' intron.info.csv > genes_list.txt
#sort -t ',' -k1,1 -V genes_list.txt > sorted_genes_list.txt

#awk 'BEGIN{FS=","} NR==FNR{genes[$1]; next} $1 ~ /^>/ {match($1, />(.*)/, m); if (m[1] in genes) {p=1} else {p=0}} p' sorted_genes_list.txt modified_headers.fasta > genes_of_interest.fasta

