#!/bin/bash
#SBATCH --job-name=azenta
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err

filename="${2%.gz}"
filename="${filename%.fastq}"
filename="${filename%.fq}"
filename="${filename%_1}"
filename="${filename%_1P}"
genome="${1%/}"
genome="${genome##*/}"


mkdir -p "${filename}_${genome}_random1_max6_result"

~/tools/STAR --runThreadN 12 \
--genomeDir "${1}" \
--readFilesIn "${2}" "${3}" \
--readFilesCommand zcat \
--outFileNamePrefix "./${filename}_${genome}_random1_max6_result/" \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 6 \
--sjdbOverhang 149 \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--limitBAMsortRAM 8000000000

cd ${filename}_${genome}_random1_max6_result/

samtools index Aligned.sortedByCoord.out.bam

/home/s1081033/miniforge3/bin/bamCoverage -bs 1 -p 12 -b Aligned.sortedByCoord.out.bam --filterRNAstrand forward -o "${filename}_${genome}_random1_max6_for.bw"
/home/s1081033/miniforge3/bin/bamCoverage -bs 1 -p 12 -b Aligned.sortedByCoord.out.bam --filterRNAstrand reverse -o "${filename}_${genome}_random1_max6_rev.bw" 
