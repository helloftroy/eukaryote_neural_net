#!/bin/bash
#SBATCH --job-name=azenta
#SBATCH --partition=cpu
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err

#usage sbatch S5_star-qc.sh $genomedir $forward.bw

genome="${1%/}"
conintronjunc="${genome}/*con.intron.junction.e1i1.bed"
allintronjunc="${genome}/*all.intron.junction.e1i1.bed"
exonbytxfile="${genome}/*exonbytx.bed"
allintronfile="${genome}/*all.intron.bed"

echo $allintronfile
#genome="${genome##*/}"

filename="${2%_for.bw}"

bigWigToBedGraph $2 ${filename}_for.bedgraph
bigWigToBedGraph ${filename}_rev.bw ${filename}_rev.bedgraph

awk -v OFS="\t" '{print $0,"ngs","+"}' ${filename}_for.bedgraph > ${filename}_for_s.bedgraph #add strand
awk -v OFS="\t" '{print $0,"ngs","-"}' ${filename}_rev.bedgraph > ${filename}_rev_s.bedgraph #add strand

cat ${filename}_for_s.bedgraph ${filename}_rev_s.bedgraph > ${filename}.bedgraph #this is the only file
rm ${filename}_for.bedgraph ${filename}_rev.bedgraph ${filename}_for_s.bedgraph ${filename}_rev_s.bedgraph

module load bedtools/2-2.29.2
bedtools intersect -a ${conintronjunc} -b ${filename}.bedgraph -wb -s | awk '{print $4,$10}' | sort  -V -s -k1,1> ${filename}.con.intron.e1i1.signal.dt
bedtools intersect -a ${allintronjunc} -b ${filename}.bedgraph -wb -s | awk '{print $4,$10}' | sort  -V -s -k1,1> ${filename}.all.intron.e1i1.signal.dt
bedtools intersect -a ${exonbytxfile} -b ${filename}.bedgraph -wb -s | awk '{print $4,$10,$3-$2}' | awk '{print $1, $2 * $3}' | awk '{a[$1] += $2} END{for (i in a) print i, a[i]}' | sort  -V -s -k1,1> ${filename}.txbyexon.signal.dt
module load R/4.2.2
Rscript  ~/tools/S5_star-qc.R  Log.final.out  ${filename}.con.intron.e1i1.signal.dt ${filename}.all.intron.e1i1.signal.dt SJ.out.tab $allintronfile
