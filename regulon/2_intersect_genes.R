library(openxlsx)

rm(list = ls())

setwd('~/project/PSVnet/V8/regulon/1_data/')
load('../../../gff.RData')

# cp -r ../chip/04_anno/chip_target .

TF = 'hrpL'
chip_name = 'hrpL'
rna_name = 'hrpL'

############################# KB ############################# 

intersect_gene <- function(TF, chip_name, rna_name) {
  print(rna_name)
  peak_occupy_gene <- tryCatch(read.csv(paste('chip_target/',chip_name,'_1.targetgenes.txt',sep = ''),header = F, sep = '\t',stringsAsFactors = F)[,1],error=function(e) NULL)

  deg_kb <- tryCatch(read.csv(paste('DEG_KB/',rna_name,'_DEGs_adjp_with_log2FC.txt',sep = ''),quote = '',header = F, sep = '\t',stringsAsFactors = F),error=function(e) NULL)
  rownames(deg_kb) <- deg_kb$V3
  
  funtional_gene_kb <- intersect(peak_occupy_gene,deg_kb$V3)
  
  funtional_gene_kb <- gff_m[gff_m$locus_tag %in% funtional_gene_kb, ]$locus_tag
  gene_name_kb <- gff_m[gff_m$locus_tag %in% funtional_gene_kb, ]$gene_name

  out <- data.frame(gene_name_kb,deg_kb[funtional_gene_kb,c(1,3)])
  
  if (length(funtional_gene_kb)>0) {
    write.table(out,file=paste('functional_target_KB/',TF,'.txt',sep = ''),row.names=F,col.names=F,sep="\t",quote=F)
  }

  c(chip_name,rna_name,length(peak_occupy_gene), nrow(deg_kb),length(funtional_gene_kb),paste(funtional_gene_kb,collapse=" "),paste(gene_name_kb,collapse=" "))
}

# TF_mutant, chip_name, rna_name

stat <- rbind(
intersect_gene('hrpL','hrpL','hrpL'),
intersect_gene('hrpR','hrpR','hrpR'),
intersect_gene('hrpS','hrpS','hrpS'),
intersect_gene('rhpS','rhpR','rhpS'),
intersect_gene('phoQ','phoP','phoQ'),
intersect_gene('mgrA','mgrA','mgrA'),
intersect_gene('tsiS','tsiR','tsiS'),
intersect_gene('pilR','pilR','pilR'),
intersect_gene('pilS','pilR','pilS'),
intersect_gene('ompR','ompR','ompR'),
intersect_gene('envZ','ompR','envZ'),
intersect_gene('gacA','gacA','gacA'),
intersect_gene('gacS','gacA','gacS'),
intersect_gene('psrA','psrA','psrA'),
intersect_gene('aefR','aefR','aefR'),
intersect_gene('vfr','vfr','vfr'),
intersect_gene('algU','algU','algU'),
intersect_gene('cvsS','cvsR','cvsS'),
intersect_gene('rpoN','rpoN','rpoN')
)

colnames(stat) <- c('TF','mutant','# peak occupied gene','# DEG KB','# intersect gene KB','intersect gene KB','gene name KB')

stat_KB <- stat

#write.xlsx(stat_KB, 'stat_KB.xlsx', col.names = T,row.names = F)
#write.table(stat_KB, 'stat_KB.csv',sep = '\t',quote = F,col.names = T,row.names = F)



############################# MM ############################# 

intersect_gene <- function(TF, chip_name, rna_name) {
  print(rna_name)
  peak_occupy_gene <- tryCatch(read.csv(paste('chip_target/',chip_name,'_1.targetgenes.txt',sep = ''),header = F, sep = '\t',stringsAsFactors = F)[,1],error=function(e) NULL)
  
  deg_mm <- tryCatch(read.csv(paste('DEG_MM/',rna_name,'_DEGs_adjp_with_log2FC.txt',sep = ''),quote = '',header = F, sep = '\t',stringsAsFactors = F),error=function(e) NULL)
  rownames(deg_mm) <- deg_mm$V3
  
  funtional_gene_mm <- intersect(peak_occupy_gene,deg_mm$V3)
  
  funtional_gene_mm <- gff_m[gff_m$locus_tag %in% funtional_gene_mm, ]$locus_tag
  gene_name_mm <- gff_m[gff_m$locus_tag %in% funtional_gene_mm, ]$gene_name
  
  out <- data.frame(gene_name_mm,deg_mm[funtional_gene_mm,c(1,3)])
  
  if (length(funtional_gene_mm)>0) {
    write.table(out,file=paste('functional_target_MM/',TF,'.txt',sep = ''),row.names=F,col.names=F,sep="\t",quote=F)
  }
  
  c(chip_name,rna_name,length(peak_occupy_gene), nrow(deg_mm),length(funtional_gene_mm),paste(funtional_gene_mm,collapse=" "),paste(gene_name_mm,collapse=" "))
}

stat <- rbind(
  intersect_gene('hrpL','hrpL','hrpL'),
  intersect_gene('hrpR','hrpR','hrpR'),
  intersect_gene('hrpS','hrpS','hrpS'),
  intersect_gene('rhpS','rhpR','rhpS'),
  intersect_gene('phoQ','phoP','phoQ'),
  intersect_gene('mgrA','mgrA','mgrA'),
  intersect_gene('tsiS','tsiR','tsiS'),
  intersect_gene('pilR','pilR','pilR'),
  intersect_gene('pilS','pilR','pilS'),
  intersect_gene('ompR','ompR','ompR'),
  intersect_gene('envZ','ompR','envZ'),
  intersect_gene('gacA','gacA','gacA'),
  intersect_gene('gacS','gacA','gacS'),
  intersect_gene('psrA','psrA','psrA'),
  intersect_gene('aefR','aefR','aefR'),
  intersect_gene('vfr','vfr','vfr'),
  intersect_gene('algU','algU','algU'),
  intersect_gene('cvsS','cvsR','cvsS'),
  intersect_gene('rpoN','rpoN','rpoN')
)

colnames(stat) <- c('TF','mutant','# peak occupied gene','# DEG MM','# intersect gene MM','intersect gene MM','gene name MM')

stat_MM = stat

#write.xlsx(stat_MM, 'stat_MM.xlsx', col.names = T,row.names = F)
#write.table(stat_MM, 'stat_MM.csv',sep = '\t',quote = F,col.names = T,row.names = F)

save(stat_KB,stat_MM, file = '2_stat.RData')
