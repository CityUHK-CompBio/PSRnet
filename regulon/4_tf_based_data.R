library(dplyr)
library(varhandle)

rm(list=ls())
setwd('~/project/PSVnet/V8/regulon/1_data/')

load('3_all_net.RData')

##############  chip-seq target ##############  
net_all <- chip_net

tf <- unique(net_all$tf_id)
stat = NULL
tf_id <- tf[1]
for (tf_id in tf){
  tf_name <- net_all[net_all$tf_id == tf_id, 'tf_name'][1]
  target_id <- net_all[net_all$tf_id == tf_id, 'gene_id']
  target_id <- unique(unfactor(target_id))
  count <- length(target_id)
  target_id <- paste(target_id, collapse = ' ')
  target_name <- net_all[net_all$tf_id == tf_id, 'gene_name']
  target_name <- unique(target_name)
  target_name <- paste(target_name, collapse = ' ')
  stat = rbind(stat, data.frame(tf_id,tf_name,count,target_id,target_name))
}
stat <- stat[with(stat, order(stat$tf_name)), ]

chip_tf <- stat


############## deg KB ##############  

# by mutant
net_all <- deg_KB_net.mutant

tf <- unique(net_all$tf_id)
stat = NULL
tf_id <- tf[1]
for (tf_id in tf){
  tf_name <- net_all[net_all$tf_id == tf_id, 'tf_name'][1]
  target_id <- net_all[net_all$tf_id == tf_id, 'gene_id']
  target_id <- unique(unfactor(target_id))
  count <- length(target_id)
  target_id <- paste(target_id, collapse = ' ')
  target_name <- net_all[net_all$tf_id == tf_id, 'gene_name']
  target_name <- unique(target_name)
  target_name <- paste(target_name, collapse = ' ')
  stat = rbind(stat, data.frame(tf_id,tf_name,count,target_id,target_name))
}
stat$tf_name <- as.vector(stat$tf_name)
stat <- stat[with(stat, order(stat$tf_name)), ]

deg_KB_tf.mutant <- stat

# by tf
net_all <- deg_KB_net

tf <- unique(net_all$tf_id)
stat = NULL
tf_id <- tf[1]
for (tf_id in tf){
  tf_name <- net_all[net_all$tf_id == tf_id, 'tf_name'][1]
  target_id <- net_all[net_all$tf_id == tf_id, 'gene_id']
  target_id <- unique(unfactor(target_id))
  count <- length(target_id)
  target_id <- paste(target_id, collapse = ' ')
  target_name <- net_all[net_all$tf_id == tf_id, 'gene_name']
  target_name <- unique(target_name)
  target_name <- paste(target_name, collapse = ' ')
  stat = rbind(stat, data.frame(tf_id,tf_name,count,target_id,target_name))
}
stat$tf_name <- as.vector(stat$tf_name)
stat <- stat[with(stat, order(stat$tf_name)), ]

deg_KB_tf <- stat


############## deg MM ##############  
# by mutant
net_all <- deg_MM_net.mutant

tf <- unique(net_all$tf_id)
stat = NULL
tf_id <- tf[1]
for (tf_id in tf){
  tf_name <- net_all[net_all$tf_id == tf_id, 'tf_name'][1]
  target_id <- net_all[net_all$tf_id == tf_id, 'gene_id']
  target_id <- unique(unfactor(target_id))
  count <- length(target_id)
  target_id <- paste(target_id, collapse = ' ')
  target_name <- net_all[net_all$tf_id == tf_id, 'gene_name']
  target_name <- unique(target_name)
  target_name <- paste(target_name, collapse = ' ')
  stat = rbind(stat, data.frame(tf_id,tf_name,count,target_id,target_name))
}
stat$tf_name <- as.vector(stat$tf_name)
stat <- stat[with(stat, order(stat$tf_name)), ]

deg_MM_tf.mutant <- stat


# by tf
net_all <- deg_MM_net

tf <- unique(net_all$tf_id)
stat = NULL
tf_id <- tf[1]
for (tf_id in tf){
  tf_name <- net_all[net_all$tf_id == tf_id, 'tf_name'][1]
  target_id <- net_all[net_all$tf_id == tf_id, 'gene_id']
  target_id <- unique(unfactor(target_id))
  count <- length(target_id)
  target_id <- paste(target_id, collapse = ' ')
  target_name <- net_all[net_all$tf_id == tf_id, 'gene_name']
  target_name <- unique(target_name)
  target_name <- paste(target_name, collapse = ' ')
  stat = rbind(stat, data.frame(tf_id,tf_name,count,target_id,target_name))
}
stat$tf_name <- as.vector(stat$tf_name)
stat <- stat[with(stat, order(stat$tf_name)), ]

deg_MM_tf <- stat



############## functional target KB ##############  
# by tf
net_all <- functional_KB_net

unique(net_all$tf_name)
tf <- unique(net_all$tf_id)
stat = NULL
tf_id <- tf[11]
for (tf_id in tf){
  tf_name <- net_all[net_all$tf_id == tf_id, 'tf_name'][1]
  target_id <- net_all[net_all$tf_id == tf_id, 'gene_id']
  target_id <- unique(unfactor(target_id))
  count <- length(target_id)
  target_id <- paste(target_id, collapse = ' ')
  target_name <- net_all[net_all$tf_id == tf_id, 'gene_name']
  target_name <- unique(target_name)
  target_name <- paste(target_name, collapse = ' ')
  stat = rbind(stat, data.frame(tf_id,tf_name,count,target_id,target_name))
}
stat$tf_name <- as.vector(stat$tf_name)
stat <- stat[with(stat, order(stat$tf_name)), ]

functional_KB_tf <- stat


############## functional target MM ##############  
# by tf
net_all <- functional_MM_net

tf <- unique(net_all$tf_id)
stat = NULL
tf_id <- tf[1]
for (tf_id in tf){
  tf_name <- net_all[net_all$tf_id == tf_id, 'tf_name'][1]
  target_id <- net_all[net_all$tf_id == tf_id, 'gene_id']
  target_id <- unique(unfactor(target_id))
  count <- length(target_id)
  target_id <- paste(target_id, collapse = ' ')
  target_name <- net_all[net_all$tf_id == tf_id, 'gene_name']
  target_name <- unique(target_name)
  target_name <- paste(target_name, collapse = ' ')
  stat = rbind(stat, data.frame(tf_id,tf_name,count,target_id,target_name))
}
stat$tf_name <- as.vector(stat$tf_name)
stat <- stat[with(stat, order(stat$tf_name)), ]

functional_MM_tf <- stat










save(chip_tf,
     deg_KB_tf.mutant,deg_KB_tf,
     deg_MM_tf.mutant,deg_MM_tf,
     functional_KB_tf,
     functional_MM_tf,
     file = '4_tf_based_data.RData')
