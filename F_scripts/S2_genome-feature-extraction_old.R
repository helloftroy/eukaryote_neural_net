#!/usr/bin/env Rscript
args = commandArgs(trailing = TRUE)

# Librarys and inputs -----------------------------------------------------
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('stringr')) install.packages('stringr'); library('stringr')
options(scipen=999)

annotation.file <- args[1]
# annotation.file <- 'corn_JHAX_v5_genome.igv.gff'

filetype <- 'gff'
if(str_detect(annotation.file,'gtf$')){
  filetype <- 'gtf'
}
filename <- str_remove_all(annotation.file,'gtf$|gff$')

### read gtf file
anno.dt <- fread(annotation.file,sep = '\t',header = F)

print(paste0(filetype,' file read at', Sys.time()))
### exon location extraction

exon.dt <- anno.dt[V3=='exon',.(V1,V4,V5,V7,V9)] %>%
  `colnames<-`(c('chr','start','end','strand','description'))

if(filetype == 'gff'){
  exon.dt[,tx_id:=str_extract(description,"Parent=.*")]
  exon.dt[,tx_id:=str_remove(tx_id,';.*')]
  exon.dt[,tx_id:=str_remove_all(tx_id,'Parent=|;')]

  tx.dt<- anno.dt[V3=='mRNA',.(V1,V4,V5,V7,V9)] %>%
    `colnames<-`(c('chr','start','end','strand','description'))

  tx.dt[,tx_id:=str_extract(description,"ID=.*")]
  tx.dt[,tx_id:=str_remove(tx_id,';.*')]
  tx.dt[,tx_id:=str_remove_all(tx_id,'ID=|;')]

  tx.dt[,gene_id:=str_extract(description,"Parent=.*")]
  tx.dt[,gene_id:=str_remove(gene_id,';.*')]
  tx.dt[,gene_id:=str_remove_all(gene_id,'Parent=|;')]

  exon.dt <- merge(exon.dt,tx.dt[,.(tx_id,gene_id)],by='tx_id')
}
if(filetype == 'gtf'){
  exon.dt[,tx_id:=str_extract(description,"transcript_id .*")]
  exon.dt[,tx_id:=str_remove(tx_id,';.*')]
  exon.dt[,tx_id:=str_remove_all(tx_id,'transcript_id |;|"')]

  exon.dt[,gene_id:=str_extract(description,"gene_id .*")]
  exon.dt[,gene_id:=str_remove(gene_id,';.*')]
  exon.dt[,gene_id:=str_remove_all(gene_id,'gene_id |;|"')]
}

exon.dt[,description := paste0(gene_id,';',tx_id)]
exon.dt[,gene_id := NULL]
exon.dt[,tx_id := NULL]

#if there is a transcript list, filter at here again
exon.bed.dt <- copy(exon.dt)
exon.bed.dt[,start:= start-1]
exon.bed.dt[,score:=0]
exon.bed.dt <- exon.bed.dt[,.(chr,start,end,description,score,strand)]
write.table(exon.bed.dt, file = paste0(filename,'exonbytx.bed'),sep = '\t',quote = F,row.names = F,col.names = F) #write table

# print(paste0(filename,'exonbytx.bed',' file wrote at', Sys.time()))

### generate intron table
intron.total.dt <- copy(exon.dt)
intron.total.dt <- intron.total.dt[order(description,start)]
# intron.dt <- intron.dt[!str_detect(chr,'B')] # for B73 annotation

#most left intron will be  first intron
intron.total.dt[,next_exon_start:=shift(start,-1), by=description]
intron.total.dt[,start:= end] #intron start equal to last exon end
intron.total.dt[,end:=next_exon_start-1] # intron end equal to next exon start -1
intron.total.dt <- intron.total.dt[!is.na(end)]
intron.total.dt[,gene_id := str_remove(description,';.*')]

intron.loc.dt <- copy(intron.total.dt)
intron.loc.dt[,description := NULL]
intron.loc.dt <- intron.loc.dt[order(gene_id,start)] %>% unique() #remove single exon and duplicated intron.

intron.loc.dt[strand == '+',intronid:=paste0(gene_id,'.intron',seq(.N)), by=gene_id] #number the introns.
intron.loc.dt[strand == '-',intronid:=paste0(gene_id,'.intron',rev(seq(.N))), by=gene_id] #number the introns.

intron.loc.dt[,score:=0]
intron.bed.dt <- intron.loc.dt[,.(chr,start,end,intronid,score,strand)]
write.table(intron.bed.dt, file = paste0(filename,'all.intron.bed'),sep = '\t',quote = F,row.names = F,col.names = F) #write table
# print(paste0(filename,'intron.bed',' file wrote at', Sys.time()))


intron.count.dt <- merge(intron.total.dt,intron.loc.dt,by=c('chr','start','end','gene_id','strand'), all=T)
intronbytx.bed.dt <- intron.count.dt[,.(chr,start,end,description,score,strand)]
intronbytx.bed.dt <- intronbytx.bed.dt[order(description,start)]
write.table(intronbytx.bed.dt, file = paste0(filename,'intronbytx.bed'),sep = '\t',quote = F,row.names = F,col.names = F) #write table

tx.count.dt <- intron.total.dt[,.(description,gene_id)]
tx.count.dt <- unique(tx.count.dt)
tx.count.dt <- tx.count.dt[,.(tx.count=.N), by=.(gene_id)]

intron.count.dt <- intron.count.dt[,.(gene_id,intronid)]
intron.count.dt <- intron.count.dt[,.(intron.count=.N), by = .(gene_id,intronid)]
intron.count.dt <- merge(intron.count.dt, tx.count.dt, by='gene_id',all = T)
intron.count.dt[,type:='constitutive']
intron.count.dt[intron.count != tx.count,type:='alternative']
table(intron.count.dt$type)
intron.consititutive.bed.dt <- intron.loc.dt[intronid %in% intron.count.dt[type == 'constitutive']$intronid]
intron.consititutive.bed.dt <- intron.consititutive.bed.dt[,.(chr,start,end,intronid,score,strand)]
write.table(intron.consititutive.bed.dt, file = paste0(filename,'con.intron.bed'),sep = '\t',quote = F,row.names = F,col.names = F) #write table
print(paste0(filename,'constitutive.intron.bed',' file wrote at', Sys.time()))

# junction table
Junction_find <- function(bedfile,exonN,intronN){
  up.intron.dt <- copy(bedfile)
  up.intron.dt[,end:=start+intronN]
  up.intron.dt[,intronid:=paste0(intronid,'.up.intron')]

  up.exon.dt <- copy(bedfile)
  up.exon.dt[,end:=start]
  up.exon.dt[,start:=start-exonN]
  up.exon.dt[,intronid:=paste0(intronid,'.up.exon')]

  down.intron.dt <- copy(bedfile)
  down.intron.dt[,start:=end - intronN]
  down.intron.dt[,intronid:=paste0(intronid,'.down.intron')]

  down.exon.dt <- copy(bedfile)
  down.exon.dt[,start:=end]
  down.exon.dt[,end:=end + exonN]
  down.exon.dt[,intronid:=paste0(intronid,'.down.exon')]

  junction.dt <- rbind(up.intron.dt,up.exon.dt,down.intron.dt,down.exon.dt)
  junction.dt[strand == '+', intronid := str_replace(intronid, 'up','sd')]
  junction.dt[strand == '+', intronid := str_replace(intronid, 'down','sa')]
  junction.dt[strand == '-', intronid := str_replace(intronid, 'up','sa')]
  junction.dt[strand == '-', intronid := str_replace(intronid, 'down','sd')]
  junction.dt <- junction.dt[order(intronid)]
  return(junction.dt)
}

con.intron.junction.e1i1.bed.dt <- Junction_find(intron.consititutive.bed.dt,1,1)
# junction.e2i2.bed.dt <- Junction_find(intron.bed.dt,2,2)
all.intron.junction.e1i1.bed.dt <- Junction_find(intron.bed.dt,1,1)
all.intron.junction.e10i12.bed.dt <- Junction_find(intron.bed.dt,10,12)

write.table(con.intron.junction.e1i1.bed.dt, file =paste0(filename,'con.intron.junction.e1i1.bed'),
            sep = '\t',quote = F,row.names = F,col.names = F) # for splicing eff calculation

write.table(all.intron.junction.e1i1.bed.dt, file =paste0(filename,'all.intron.junction.e1i1.bed'),
            sep = '\t',quote = F,row.names = F,col.names = F) # for splicing feature annotation

write.table(all.intron.junction.e10i12.bed.dt, file =paste0(filename,'all.intron.junction.e10i12.bed'),
            sep = '\t',quote = F,row.names = F,col.names = F) # for splicing feature annotation



