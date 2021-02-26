
library(dplyr)
library(varhandle)
library(openxlsx)

rm(list=ls())
setwd('~/project/PSVnet/V8/regulon/1_data/')

load('3_all_net.RData')

##############  chip-seq target ##############  
net_all <- chip_net

target <- unique(unfactor(net_all$gene_id))
stat = NULL
target_id <- target[1]
for (target_id in target){
  target_name <- net_all[net_all$gene_id == target_id, 'gene_name'][1]
  tf_id <- net_all[net_all$gene_id == target_id, 'tf_id']
  tf_id <- unique(unfactor(tf_id))
  count <- length(tf_id)
  tf_id <- paste(tf_id, collapse = ' ')
  tf_name <- net_all[net_all$gene_id == target_id, 'tf_name']
  tf_name <- unique(unfactor(tf_name))
  tf_name <- paste(tf_name, collapse = ' ')
  stat = rbind(stat, data.frame(target_id,target_name,count,tf_id,tf_name))
}

chip_target <- stat
write.xlsx(chip_target[chip_target$count > 1,], 'chip_target.xlsx')


############## deg KB ##############  

# by mutant
net_all <- deg_KB_net.mutant

target <- unique(unfactor(net_all$gene_id))
stat = NULL
target_id <- target[1]
for (target_id in target){
  target_name <- net_all[net_all$gene_id == target_id, 'gene_name'][1]
  tf_id <- net_all[net_all$gene_id == target_id, 'tf_id']
  tf_id <- unique(tf_id)
  count <- length(tf_id)
  tf_id <- paste(tf_id, collapse = ' ')
  tf_name <- net_all[net_all$gene_id == target_id, 'tf_name']
  tf_name <- unique(tf_name)
  tf_name <- paste(tf_name, collapse = ' ')
  logFC = net_all[net_all$gene_id == target_id, 'logFC']
  logFC = logFC[which.max( abs(logFC) )]
  stat = rbind(stat, data.frame(target_id,target_name,logFC,count,tf_id,tf_name))
}

deg_KB_target.mutant <- stat



# by tf
net_all <- deg_KB_net

target <- unique(unfactor(net_all$gene_id))
stat = NULL
target_id <- target[1]
for (target_id in target){
  target_name <- net_all[net_all$gene_id == target_id, 'gene_name'][1]
  tf_id <- net_all[net_all$gene_id == target_id, 'tf_id']
  tf_id <- unique(tf_id)
  count <- length(tf_id)
  tf_id <- paste(tf_id, collapse = ' ')
  tf_name <- net_all[net_all$gene_id == target_id, 'tf_name']
  tf_name <- unique(tf_name)
  tf_name <- paste(tf_name, collapse = ' ')
  logFC = net_all[net_all$gene_id == target_id, 'logFC']
  logFC = logFC[which.max( abs(logFC) )]
  stat = rbind(stat, data.frame(target_id,target_name,logFC,count,tf_id,tf_name))
}

deg_KB_target <- stat
write.xlsx(deg_KB_target[deg_KB_target$count > 1,], 'deg_KB_target.xlsx')



############## deg MM ##############  

# by mutant
net_all <- deg_MM_net.mutant

target <- unique(unfactor(net_all$gene_id))
stat = NULL
target_id <- target[1]
for (target_id in target){
  target_name <- net_all[net_all$gene_id == target_id, 'gene_name'][1]
  tf_id <- net_all[net_all$gene_id == target_id, 'tf_id']
  tf_id <- unique(tf_id)
  count <- length(tf_id)
  tf_id <- paste(tf_id, collapse = ' ')
  tf_name <- net_all[net_all$gene_id == target_id, 'tf_name']
  tf_name <- unique(tf_name)
  tf_name <- paste(tf_name, collapse = ' ')
  logFC = net_all[net_all$gene_id == target_id, 'logFC']
  logFC = logFC[which.max( abs(logFC) )]
  stat = rbind(stat, data.frame(target_id,target_name,logFC,count,tf_id,tf_name))
}

deg_MM_target.mutant <- stat


# by tf
net_all <- deg_MM_net

target <- unique(unfactor(net_all$gene_id))
stat = NULL
target_id <- target[1]
for (target_id in target){
  target_name <- net_all[net_all$gene_id == target_id, 'gene_name'][1]
  tf_id <- net_all[net_all$gene_id == target_id, 'tf_id']
  tf_id <- unique(tf_id)
  count <- length(tf_id)
  tf_id <- paste(tf_id, collapse = ' ')
  tf_name <- net_all[net_all$gene_id == target_id, 'tf_name']
  tf_name <- unique(tf_name)
  tf_name <- paste(tf_name, collapse = ' ')
  logFC = net_all[net_all$gene_id == target_id, 'logFC']
  logFC = logFC[which.max( abs(logFC) )]
  stat = rbind(stat, data.frame(target_id,target_name,logFC,count,tf_id,tf_name))
}

deg_MM_target <- stat
write.xlsx(deg_MM_target[deg_MM_target$count > 1,], 'deg_MM_target.xlsx')


############## functional target KB ##############  

# by tf
net_all <- functional_KB_net

target <- unique(unfactor(net_all$gene_id))
stat = NULL
target_id <- target[1]
for (target_id in target){
  target_name <- net_all[net_all$gene_id == target_id, 'gene_name'][1]
  tf_id <- net_all[net_all$gene_id == target_id, 'tf_id']
  tf_id <- unique(tf_id)
  count <- length(tf_id)
  tf_id <- paste(tf_id, collapse = ' ')
  tf_name <- net_all[net_all$gene_id == target_id, 'tf_name']
  tf_name <- unique(tf_name)
  tf_name <- paste(tf_name, collapse = ' ')
  direct = net_all[net_all$gene_id == target_id, 'direct']
  direct = direct[which.max( abs(direct) )]
  stat = rbind(stat, data.frame(target_id,target_name,direct,count,tf_id,tf_name))
}

functional_KB_target <- stat
write.xlsx(functional_KB_target[functional_KB_target$count > 1,], 'functional_KB_target.xlsx')



############## functional target MM ##############  

# by tf
net_all <- functional_MM_net

target <- unique(unfactor(net_all$gene_id))
stat = NULL
target_id <- target[1]
for (target_id in target){
  target_name <- net_all[net_all$gene_id == target_id, 'gene_name'][1]
  tf_id <- net_all[net_all$gene_id == target_id, 'tf_id']
  tf_id <- unique(tf_id)
  count <- length(tf_id)
  tf_id <- paste(tf_id, collapse = ' ')
  tf_name <- net_all[net_all$gene_id == target_id, 'tf_name']
  tf_name <- unique(tf_name)
  tf_name <- paste(tf_name, collapse = ' ')
  direct = net_all[net_all$gene_id == target_id, 'direct']
  direct = direct[which.max( abs(direct) )]
  stat = rbind(stat, data.frame(target_id,target_name,direct,count,tf_id,tf_name))
}

functional_MM_target <- stat
write.xlsx(functional_MM_target[functional_MM_target$count > 1,], 'functional_MM_target.xlsx')



save(chip_target,
     deg_KB_target.mutant,deg_KB_target,
     deg_MM_target.mutant,deg_MM_target,
     functional_KB_target,
     functional_MM_target,
     file = '5_target_based_data.RData')



