#!/bin/bash
#SBATCH --job-name=starqc
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err

#usage qsub S6_intron_info.sh "*QCfiles" ".*junction.csv.files" ".*splicing.csv.files" $outputname
#note that $2-3 need use R version wild card such as .*
outputname=$4

for i in `ls $1`
do
        filename=${i%%_genome*}
        sed '/^$/d' $i | awk -F "|" '{print $1}' | tr '\n' ',' | awk -v FS="," -v OFS="," '{print "filename",$0}'>> ${outputname}.star.qc.tb
        sed '/^$/d' $i | awk -F "|" '{print $2}' | tr '\n' ',' | awk -v FS="," -v OFS=","  -v name=$filename '{print name,$0}' >> ${outputname}.star.qc.tb
done
sed -n '1p;n;p' ${outputname}.star.qc.tb | sed 's/,$//g' > ${outputname}.star.qc.table.csv
rm ${outputname}.star.qc.tb

module load R/4.2.2
Rscript  ~/tools/S6_intron_info.R $2  $3  $4

