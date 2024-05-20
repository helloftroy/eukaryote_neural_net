#!/bin/bash -l
set  -euo pipefail
#$ -cwd
#$ -j y
#$ -N 'trim-pe'
#$ -pe smp 4
#$ -l m_mem_free=9G

filename="${1%.gz}"
filename="${filename%.fastq}"
filename="${filename%.fq}"
filename="${filename%_1}"

cp -u ~/genomes/com.adaptors.fa ./

java -jar ~/tools/trimmomatic.jar PE -threads 4 \
-basein ${1} \
-baseout "${filename}.fq.gz" \
ILLUMINACLIP:com.adaptors.fa:2:30:10:3:true LEADING:26 TRAILING:26 MINLEN:25

rm ${filename}_1U.fq.gz ${filename}_2U.fq.gz

