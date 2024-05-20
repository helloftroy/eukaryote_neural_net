#!/usr/bin/env Rscript
args = commandArgs(trailing = TRUE)

# Librarys and inputs -----------------------------------------------------
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('stringr')) install.packages('stringr'); library('stringr')
options(scipen=999)

junction.files.pattern <- args[1]
splicing.files.pattern <- args[2]
output.filename <- args[3]

# junction.files.pattern <- "*intron.junc.csv"
junction.files <- list.files(pattern = junction.files.pattern)
# 
# splicing.files.pattern <- "*intron.splicing.csv"
splicing.files <- list.files(pattern = splicing.files.pattern)
# output.filename <- 'test'
print(junction.files)
print(splicing.files)
# read junction files ------------------------------------------------------

junction.file.list <- list()

for(i in 1:length(junction.files)){
  junction.file.list[[i]] <- fread(junction.files[i])
}

junction.dt <- rbindlist(junction.file.list)

junction.dt <- junction.dt[,.(total.count.exon.sd = sum(count.exon.sd),
                              total.count.intron.sd = sum(count.intron.sd),
                              total.count.exon.sa = sum(count.exon.sa),
                              total.count.intron.sa = sum(count.intron.sa)),
                           by=.(intron_id)]


junction.dt[total.count.exon.sd > 4,sd.eff:=round(100*(total.count.exon.sd - total.count.intron.sd)/total.count.exon.sd)]
junction.dt[total.count.exon.sa > 4,sa.eff:=round(100*(total.count.exon.sa - total.count.intron.sa)/total.count.exon.sa)]

junction.dt[sd.eff<0, sd.eff:=0]
junction.dt[sa.eff<0, sa.eff:=0]
#write.csv(junction.dt, file=paste0(output.filename,'.combine.intron.efficiency.csv'), row.names = F, quote = F)
print('ok')
# read splicing files ------------------------------------------------------

splicing.file.list <- list()

for(i in 1:length(splicing.files)){
  splicing.file.list[[i]] <- fread(splicing.files[i])
}

splicing.dt <- rbindlist(splicing.file.list)

splicing.dt <- splicing.dt[,.(total.correct.splicing.reads = sum(correct.splicing.reads),
                              total.mis.3ss.reads = sum(mis.3ss.reads),
                              total.mis.5ss.reads = sum(mis.5ss.reads)),
                           by=.(intron.id)]
splicing.dt[,sum.reads := total.correct.splicing.reads + total.mis.3ss.reads + total.mis.5ss.reads]
splicing.dt[,accuracy := round(total.correct.splicing.reads/sum.reads *100)]

#write.csv(splicing.dt, file=paste0(output.filename,'.combine.intron.accuracy.csv'), row.names = F, quote = F)

# write table -------------------------------------------------------------
intron.info.dt <- merge(junction.dt,splicing.dt,by.x='intron_id',by.y='intron.id',all=T)

write.csv(intron.info.dt, file=paste0(output.filename,'.combine.intron.info.csv'), row.names = F, quote = F)

