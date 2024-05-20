#!/bin/bash -l
set  -euo pipefail
#$ -now y
#$ -clear
#$ -cwd
#$ -j y
#$ -N 'feature'
#$ -pe smp 1
#$ -l m_mem_free=5G

#usage: qsub S2_genome-feature-extraction.sh $genome.fasta $genome.gtf/gff

module load R/4.2.2

filename=${2%.gtf}
filename=${filename%.gff}

sed -i '/^#/d' $2
Rscript ~/tools/S2_genome-feature-extraction.R $2

module load bedtools/2-2.29.2
bedtools getfasta -s -fi $1 -bed ${filename}.all.intron.bed > ${filename}.all.intron.fasta
bedtools getfasta -s -fi $1 -bed ${filename}.all.intron.junction.e10i12.bed | sed '/^>/d'> ${filename}.all.intron.junction.e10i12.seq

e10i12bed=${filename}.all.intron.junction.e10i12.bed
e10i12seq=${filename}.all.intron.junction.e10i12.seq

#-----------
cat $e10i12bed | awk '{print $4}' | sed 's/.[a-z]\+$//' | uniq | sed 's/^/>/' > junc.id.temp
sed -n 'p;n' junc.id.temp > junc.sa.id.temp
sed -n 'n;p' junc.id.temp > junc.sd.id.temp

cat $e10i12seq | sed -n 'p;n;n;n' > sa.exon.temp
cat $e10i12seq | sed -n 'n;p;n;n' > sa.intron.temp
cat $e10i12seq | sed -n 'n;n;p;n' > sd.exon.temp
cat $e10i12seq | sed -n 'n;n;n;p' > sd.intron.temp

paste -d "" sd.exon.temp sd.intron.temp > sd.seq.temp
paste -d "" sa.intron.temp sa.exon.temp > sa.seq.temp

paste -d "\\n" junc.sa.id.temp sa.seq.temp > junc.sa.fasta.temp
paste -d "\\n" junc.sd.id.temp sd.seq.temp > junc.sd.fasta.temp

echo 'ok'
#SD
module load viennaRNA/2.5.1
RNAfold -p --noPS --noDP -T 28 --MEA junc.sd.fasta.temp > junc.sd.MEA.temp
sed -n '1~7p' junc.sd.MEA.temp | sed -r 's/>//g' > junc.sd.name.temp
sed -n '3~7p' junc.sd.MEA.temp | sed 's/ .*$//' > junc.sd.MFE.struc.temp
sed -n '3~7p' junc.sd.MEA.temp | sed 's/^.* //;s/(//;s/)//' > junc.sd.MFE.score.temp
sed -n '5~7p' junc.sd.MEA.temp | sed 's/ .*$//' > junc.sd.centroid.struc.temp

paste -d ',' junc.sd.name.temp junc.sd.MFE.struc.temp junc.sd.centroid.struc.temp junc.sd.MFE.score.temp  > junc.sd.csv.temp
sed -i '1i sd.id,sd.MFEstructure,sd.Centroid.struc,sd.MFEscore' junc.sd.csv.temp

#SA
RNAfold -p --noPS --noDP -T 28 --MEA junc.sa.fasta.temp > junc.sa.MEA.temp
sed -n '1~7p' junc.sa.MEA.temp | sed -r 's/>//g' > junc.sa.name.temp
sed -n '3~7p' junc.sa.MEA.temp | sed 's/ .*$//' > junc.sa.MFE.struc.temp
sed -n '3~7p' junc.sa.MEA.temp | sed 's/^.* //;s/(//;s/)//' > junc.sa.MFE.score.temp
sed -n '5~7p' junc.sa.MEA.temp | sed 's/ .*$//' > junc.sa.centroid.struc.temp

paste -d ',' junc.sa.name.temp junc.sa.MFE.struc.temp junc.sa.centroid.struc.temp junc.sa.MFE.score.temp  > junc.sa.csv.temp
sed -i '1i sa.id,sa.MFEstructure,sa.Centroid.struc,sa.MFEscore' junc.sa.csv.temp

paste -d "," junc.sd.csv.temp junc.sa.csv.temp > ${filename}.all.intron.junction.e10i12.csv

###now generate a control which is shifting all junctions to upstream 22bp. so the sd upstream will be in exon and sa upstream will be in intron. I want to see whether there are any difference in structure.

bedtools shift -i ${filename}.all.intron.junction.e10i12.bed -g chrNameLength.txt -p -22 -m 22 > ${filename}.all.intron.junction.e10i12up22.bed
bedtools getfasta -s -fi $1 -bed ${filename}.all.intron.junction.e10i12up22.bed | sed '/^>/d'> ${filename}.all.intron.junction.e10i12up22.seq

e10i12up22bed=${filename}.all.intron.junction.e10i12up22.bed
e10i12up22seq=${filename}.all.intron.junction.e10i12up22.seq

cat $e10i12up22bed | awk '{print $4}' | sed 's/.[a-z]\+$//' | uniq | sed 's/^/>/' > junc.id.temp
sed -n 'p;n' junc.id.temp > junc.sa.id.temp
sed -n 'n;p' junc.id.temp > junc.sd.id.temp

cat $e10i12up22seq | sed -n 'p;n;n;n' > sa.exon.temp
cat $e10i12up22seq | sed -n 'n;p;n;n' > sa.intron.temp
cat $e10i12up22seq | sed -n 'n;n;p;n' > sd.exon.temp
cat $e10i12up22seq | sed -n 'n;n;n;p' > sd.intron.temp

paste -d "" sd.exon.temp sd.intron.temp > sd.seq.temp
paste -d "" sa.intron.temp sa.exon.temp > sa.seq.temp

paste -d "\\n" junc.sa.id.temp sa.seq.temp > junc.sa.fasta.temp
paste -d "\\n" junc.sd.id.temp sd.seq.temp > junc.sd.fasta.temp

echo 'ok'
#SD
RNAfold -p --noPS --noDP -T 28 --MEA junc.sd.fasta.temp > junc.sd.MEA.temp
sed -n '1~7p' junc.sd.MEA.temp | sed -r 's/>//g' > junc.sd.name.temp
sed -n '3~7p' junc.sd.MEA.temp | sed 's/ .*$//' > junc.sd.MFE.struc.temp
sed -n '3~7p' junc.sd.MEA.temp | sed 's/^.* //;s/(//;s/)//' > junc.sd.MFE.score.temp
sed -n '5~7p' junc.sd.MEA.temp | sed 's/ .*$//' > junc.sd.centroid.struc.temp

paste -d ',' junc.sd.name.temp junc.sd.MFE.struc.temp junc.sd.centroid.struc.temp junc.sd.MFE.score.temp  > junc.sd.csv.temp
sed -i '1i sd.id,sd.MFEstructure,sd.Centroid.struc,sd.MFEscore' junc.sd.csv.temp

#SA
RNAfold -p --noPS --noDP -T 28 --MEA junc.sa.fasta.temp > junc.sa.MEA.temp
sed -n '1~7p' junc.sa.MEA.temp | sed -r 's/>//g' > junc.sa.name.temp
sed -n '3~7p' junc.sa.MEA.temp | sed 's/ .*$//' > junc.sa.MFE.struc.temp
sed -n '3~7p' junc.sa.MEA.temp | sed 's/^.* //;s/(//;s/)//' > junc.sa.MFE.score.temp
sed -n '5~7p' junc.sa.MEA.temp | sed 's/ .*$//' > junc.sa.centroid.struc.temp

paste -d ',' junc.sa.name.temp junc.sa.MFE.struc.temp junc.sa.centroid.struc.temp junc.sa.MFE.score.temp  > junc.sa.csv.temp
sed -i '1i sa.id,sa.MFEstructure,sa.Centroid.struc,sa.MFEscore' junc.sa.csv.temp

paste -d "," junc.sd.csv.temp junc.sa.csv.temp > ${filename}.all.intron.junction.e10i12up22.csv


rm *temp
rm ${filename}.all.intron.junction.e10i12up22.seq
rm ${filename}.all.intron.junction.e10i12up22.bed

