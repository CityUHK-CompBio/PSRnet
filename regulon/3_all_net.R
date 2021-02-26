## creat 20199019
## get all connections from data files
## chipseq target, DEG MM, DEG KB from data
## functional target KB, functional target MM, reverse directions, add known information 

## rpon direction?


library(dplyr)
library(varhandle)
library(readxl)


rm(list=ls())
setwd('~/project/PSVnet/V8/regulon/1_data/')

load('~/project/PSVnet/gff.RData')
tf_file <- read.csv('~/project/PSVnet/tf_list.txt',sep = '\t',header = T)
mutant_file <- read.csv('~/project/PSVnet/mutant_list.txt',sep = '\t',header = T)


##############  chip-seq target ##############  

flist <- Sys.glob("chip_target/*_1.targetgenes.txt")
i = 1
net = NULL
for (i in 1 : length(flist)){
  tf = strsplit(flist[i],split = '_1')[[1]]
  tf = strsplit(tf,split = '/',2)[[1]][2]
  print(tf)
  net[[i]] <- read.table(flist[i],header=F)
  net[[i]][,'tf'] <- tf
}
net_all <- do.call("rbind",net)
net_all <- distinct(net_all)

net_all <- merge(net_all, tf_file[,c(1:3)],by.x = 'tf', by.y = 'tf_name2')
net_all <- merge(net_all, gff_m[,c(11,12)],by.x = 'V1', by.y = 'locus_tag')
net_all <- net_all[c('tf_id','tf_name','V1','gene_name')]
colnames(net_all) <- c('tf_id','tf_name','gene_id','gene_name')

chip_net <- net_all


############## deg KB ##############  

flist <- Sys.glob("DEG_KB/*")
net <- NULL
i = 1
for (i in 1 : length(flist)){
  tf = strsplit(flist[i],split = '/')[[1]][2]
  tf = strsplit(tf,split = '_')[[1]][1]
  print(tf)
  net[[i]] <- read.table(flist[i],header=F)
  net[[i]][,'tf'] <- tf
}
net_all <- do.call("rbind",net)
net_all <- distinct(net_all)

net_all <- merge(net_all, gff_m[,c(11:13)],by.x = 'tf', by.y = 'gene_name2')
net_all <- merge(net_all, gff_m[,c(11,12)],by.x = 'V3', by.y = 'locus_tag')
net_all <- net_all[c('locus_tag','gene_name.x','V3','gene_name.y','V1')]
colnames(net_all) <- c('tf_id','tf_name','gene_id','gene_name','logFC')

net_all$direct <- -net_all$logFC
net_all[net_all$tf_name == 'RpoN','direct'] = -net_all[net_all$tf_name == 'RpoN','direct'] 

deg_KB_net.mutant <- net_all

i = 12
for (i in 1:nrow(mutant_file)) {
  net_all$tf_name[net_all$tf_name == as.character(mutant_file[i,'mutant_name']) ] <- as.character(mutant_file[i,'tf_name'])
  net_all$tf_id[net_all$tf_id == as.character(mutant_file[i,'mutant_id']) ] <- as.character(mutant_file[i,'tf_id'])
}
deg_KB_net <- net_all

unique(net_all$tf_name)

############## deg MM ##############  

flist <- Sys.glob("DEG_MM/*")
net <- NULL
i = 1
for (i in 1 : length(flist)){
  tf = strsplit(flist[i],split = '/')[[1]][2]
  tf = strsplit(tf,split = '_')[[1]][1]
  print(tf)
  net[[i]] <- read.table(flist[i],header=F)
  net[[i]][,'tf'] <- tf
}
net_all <- do.call("rbind",net)
net_all <- distinct(net_all)

net_all <- merge(net_all, gff_m[,c(11:13)],by.x = 'tf', by.y = 'gene_name2')
net_all <- merge(net_all, gff_m[,c(11,12)],by.x = 'V3', by.y = 'locus_tag')
net_all <- net_all[c('locus_tag','gene_name.x','V3','gene_name.y','V1')]
colnames(net_all) <- c('tf_id','tf_name','gene_id','gene_name','logFC')

net_all$direct <- -net_all$logFC
net_all[net_all$tf_name == 'RpoN','direct'] = -net_all[net_all$tf_name == 'RpoN','direct'] 

deg_MM_net.mutant <- net_all

for (i in 1:nrow(mutant_file)) {
  net_all$tf_name[net_all$tf_name == as.character(mutant_file[i,'mutant_name']) ] <- as.character(mutant_file[i,'tf_name'])
  net_all$tf_id[net_all$tf_id == as.character(mutant_file[i,'mutant_id']) ] <- as.character(mutant_file[i,'tf_id'])
}
deg_MM_net <- net_all

#temp = known$target_name

############## known functional targets ############## 

known <- read_excel('~/project/PSVnet/known_connection.xlsx',col_names = F)
known <- known[,c(1:4)]
colnames(known) <- c('type','tf_name','target_name','direct')
known <- merge(known,gff_m[,c("gene_name","locus_tag")],by.x = 'target_name',by.y = 'gene_name')
known <- merge(known,gff_m[,c("gene_name","locus_tag")],by.x = 'tf_name',by.y = 'gene_name')
known = known[c('type',"locus_tag.y",'tf_name',"locus_tag.x",'target_name','direct')]
colnames(known) <- c('type','tf_id','tf_name','gene_id','gene_name','direct')
known$logFC <- 0
known$direct[known$direct == 1] <- 10
known$direct[known$direct == -1] <- -10

known_KB <- known[known$type == 'KB',]
known_KB <- known_KB[,c('tf_id','tf_name','gene_id','gene_name',"logFC",'direct')]
known_MM <- known[known$type == 'MM',]
known_MM <- known_MM[,c('tf_id','tf_name','gene_id','gene_name',"logFC",'direct')]


############## functional target KB ##############  

flist <- Sys.glob("functional_target_KB/*.txt")
net <- NULL
i = 1
for (i in 1 : length(flist)){
  tf = strsplit(flist[i],split = '/')[[1]][2]
  tf = strsplit(tf,split = '.txt')[[1]][1]
  print(tf)
  net[[i]] <- read.table(flist[i],header=F)
  net[[i]][,'tf'] <- tf
}
net_all <- do.call("rbind",net)
net_all <- distinct(net_all)

net_all <- merge(net_all, gff_m[,c(11:13)],by.x = 'tf', by.y = 'gene_name2')
net_all <- merge(net_all, gff_m[,c(11,12)],by.x = 'V3', by.y = 'locus_tag')
net_all <- net_all[c('locus_tag','gene_name.x','V3','gene_name.y','V2')]
colnames(net_all) <- c('tf_id','tf_name','gene_id','gene_name','logFC')

net_all$direct <- -net_all$logFC
net_all[net_all$tf_name == 'RpoN','direct'] = -net_all[net_all$tf_name == 'RpoN','direct'] 

#functional_KB_net.mutant <- net_all

for (i in 1:nrow(mutant_file)) {
  net_all$tf_name[net_all$tf_name == as.character(mutant_file[i,'mutant_name']) ] <- as.character(mutant_file[i,'tf_name'])
  net_all$tf_id[net_all$tf_id == as.character(mutant_file[i,'mutant_id']) ] <- as.character(mutant_file[i,'tf_id'])
}

functional_KB_net.withoutknown <- net_all

net_all <- rbind(net_all,known_KB)
net_all <- distinct(net_all)

self <- net_all[which(net_all$tf_name == net_all$gene_name),]
self[self$tf_name == 'RhpR', 'direct'] = 10
self <- self[self$direct == 10 | self$direct == -10,]

temp <- net_all[which(net_all$tf_name != net_all$gene_name),]

net_all <- rbind(self, temp)

functional_KB_net <- net_all


############## functional target MM ##############  

flist <- Sys.glob("functional_target_MM/*.txt")
net <- NULL
i = 1
for (i in 1 : length(flist)){
  tf = strsplit(flist[i],split = '/')[[1]][2]
  tf = strsplit(tf,split = '.txt')[[1]][1]
  print(tf)
  net[[i]] <- read.table(flist[i],header=F)
  net[[i]][,'tf'] <- tf
}
net_all <- do.call("rbind",net)
net_all <- distinct(net_all)

net_all <- merge(net_all, gff_m[,c(11:13)],by.x = 'tf', by.y = 'gene_name2')
net_all <- merge(net_all, gff_m[,c(11,12)],by.x = 'V3', by.y = 'locus_tag')
net_all <- net_all[c('locus_tag','gene_name.x','V3','gene_name.y','V2')]
colnames(net_all) <- c('tf_id','tf_name','gene_id','gene_name','logFC')

net_all$direct <- -net_all$logFC
net_all[net_all$tf_name == 'RpoN','direct'] = -net_all[net_all$tf_name == 'RpoN','direct'] 

functional_MM_net.mutant <- net_all

for (i in 1:nrow(mutant_file)) {
  net_all$tf_name[net_all$tf_name == as.character(mutant_file[i,'mutant_name']) ] <- as.character(mutant_file[i,'tf_name'])
  net_all$tf_id[net_all$tf_id == as.character(mutant_file[i,'mutant_id']) ] <- as.character(mutant_file[i,'tf_id'])
}

functional_MM_net.withoutknown <- net_all

net_all <- rbind(net_all,known_MM)
net_all <- distinct(net_all)

self <- net_all[which(net_all$tf_name == net_all$gene_name),]
#self[self$tf_name == 'RhpR', 'direct'] = 10
self[self$tf_name == 'HrpL', 'direct'] = -10
self[self$tf_name == 'PsrA', 'direct'] = -10
self <- self[self$direct == 10 | self$direct == -10,]

temp <- net_all[which(net_all$tf_name != net_all$gene_name),]

net_all <- rbind(self, temp)

functional_MM_net <- net_all




###############

save(chip_net,
     deg_KB_net.mutant,deg_KB_net,
     deg_MM_net.mutant,deg_MM_net,
     functional_KB_net,
     functional_MM_net,
     file = '3_all_net.RData')

save(known,
     functional_KB_net,functional_KB_net.withoutknown,known_KB,
     functional_MM_net,functional_MM_net.withoutknown,known_MM,
     file = '3_functional_net.RData')


curve_arrow <- functional_KB_net[which(functional_KB_net$tf_name == functional_KB_net$gene_name),]
View(curve_arrow)
curve_arrow <- functional_MM_net[which(functional_MM_net$tf_name == functional_MM_net$gene_name),]
View(curve_arrow)

