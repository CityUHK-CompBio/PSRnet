library(dplyr)
library(varhandle)

rm(list=ls())
setwd('~/project/PSVnet/V8/regulon/1_data/')

load('3_all_net.RData')




############## intersect ##############  
############## functional target KB ##############  
# by tf
net_all <- functional_KB_net

unique(net_all$tf_name)
tf <- sort(unique(net_all$tf_name))
target_list = list()

stat = NULL
tf_name <- tf[11]
for (tf_name in tf){
  tf_name <- net_all[net_all$tf_name == tf_name, 'tf_name'][1]
  target_id <- net_all[net_all$tf_name == tf_name, 'gene_id']
  target_id <- unique(unfactor(target_id))
  target_list[[tf_name]] = target_id
}

combination <- as.data.frame(t(combn(tf,2)))

i = 58
for (i in 1:nrow(combination)) {
  target_1 = target_list[[as.character(combination[i,1])]]
  target_2 = target_list[[as.character(combination[i,2])]]
  combination[i,3] = length(intersect(target_1,target_2))
  combination[i,4] = as.character(paste(intersect(target_1,target_2),collapse = ' '))
}
combination$V5 = log2(combination$V3)

functional_KB_tf.intersect <- combination


## select by frequence
select_link = NULL

tf_name <- tf[1]
for (tf_name in tf){
  df_sub = combination[(combination$V1 == tf_name)|(combination$V2 == tf_name) , ]
  df_sub = df_sub[order(df_sub$V3,decreasing = T),]
  select_link = rbind(select_link, df_sub[1,])
}
select_link$width = log2(select_link$V3)

functional_KB_tf.select_link = select_link



############## intersect ##############  
############## functional target MM ##############  
# by tf
net_all <- functional_MM_net

unique(net_all$tf_name)
tf <- unique(net_all$tf_name)
target_list = list()

stat = NULL
tf_name <- tf[11]
for (tf_name in tf){
  tf_name <- net_all[net_all$tf_name == tf_name, 'tf_name'][1]
  target_id <- net_all[net_all$tf_name == tf_name, 'gene_id']
  target_id <- unique(unfactor(target_id))
  target_list[[tf_name]] = target_id
}

combination <- as.data.frame(t(combn(tf,2)))

i = 1
for (i in 1:nrow(combination)) {
  target_1 = target_list[[as.character(combination[i,1])]]
  target_2 = target_list[[as.character(combination[i,2])]]
  combination[i,3] = length(intersect(target_1,target_2))
  combination[i,4] = as.character(paste(intersect(target_1,target_2),collapse = ' '))
}
combination$V5 = log2(combination$V3)

functional_MM_tf.intersect <- combination


## select by frequence
select_link = NULL

tf_name <- tf[1]
for (tf_name in tf){
  df_sub = combination[(combination$V1 == tf_name)|(combination$V2 == tf_name) , ]
  df_sub = df_sub[order(df_sub$V3,decreasing = T),]
  select_link = rbind(select_link, df_sub[1,])
}
select_link$width = log2(select_link$V3)

functional_MM_tf.select_link = select_link





############## intersect ##############  
############## deg KB ##############  
# by tf
net_all <- deg_KB_net

unique(net_all$tf_name)
tf <- sort(unique(net_all$tf_name))
target_list = list()

stat = NULL
tf_name <- tf[11]
for (tf_name in tf){
  tf_name <- net_all[net_all$tf_name == tf_name, 'tf_name'][1]
  target_id <- net_all[net_all$tf_name == tf_name, 'gene_id']
  target_id <- unique(unfactor(target_id))
  target_list[[tf_name]] = target_id
}

combination <- as.data.frame(t(combn(tf,2)))

i = 58
for (i in 1:nrow(combination)) {
  target_1 = target_list[[as.character(combination[i,1])]]
  target_2 = target_list[[as.character(combination[i,2])]]
  combination[i,3] = length(intersect(target_1,target_2))
  combination[i,4] = as.character(paste(intersect(target_1,target_2),collapse = ' '))
}
#combination$V4 = log2(combination$V3)

deg_KB_tf.intersect <- combination




############## intersect ##############  
############## deg MM ##############  
# by tf
net_all <- deg_MM_net

unique(net_all$tf_name)
tf <- sort(unique(net_all$tf_name))
target_list = list()

stat = NULL
tf_name <- tf[11]
for (tf_name in tf){
  tf_name <- net_all[net_all$tf_name == tf_name, 'tf_name'][1]
  target_id <- net_all[net_all$tf_name == tf_name, 'gene_id']
  target_id <- unique(unfactor(target_id))
  target_list[[tf_name]] = target_id
}

combination <- as.data.frame(t(combn(tf,2)))

i = 58
for (i in 1:nrow(combination)) {
  target_1 = target_list[[as.character(combination[i,1])]]
  target_2 = target_list[[as.character(combination[i,2])]]
  combination[i,3] = length(intersect(target_1,target_2))
  combination[i,4] = as.character(paste(intersect(target_1,target_2),collapse = ' '))
}
#combination$V4 = log2(combination$V3)

deg_MM_tf.intersect <- combination




############## intersect ##############  
############## chip ##############  
# by tf
net_all <- chip_net

unique(net_all$tf_name)
tf <- unfactor(sort(unique(net_all$tf_name)))
target_list = list()

stat = NULL
tf_name <- tf[11]
for (tf_name in tf){
  tf_name <- as.character(net_all[net_all$tf_name == tf_name, 'tf_name'][1])
  target_id <- net_all[net_all$tf_name == tf_name, 'gene_id']
  target_id <- unique(unfactor(target_id))
  target_list[[tf_name]] = target_id
}

combination <- as.data.frame(t(combn(tf,2)))

i = 58
for (i in 1:nrow(combination)) {
  target_1 = target_list[[as.character(combination[i,1])]]
  target_2 = target_list[[as.character(combination[i,2])]]
  combination[i,3] = length(intersect(target_1,target_2))
  combination[i,4] = as.character(paste(intersect(target_1,target_2),collapse = ' '))
}
#combination$V4 = log2(combination$V3)

chip_tf.intersect <- combination




save(functional_KB_tf.intersect,functional_MM_tf.intersect,
     functional_KB_tf.select_link,functional_MM_tf.select_link,
     deg_KB_tf.intersect,deg_MM_tf.intersect,
     chip_tf.intersect,
     file = '7_tf_target_intersect.RData')










