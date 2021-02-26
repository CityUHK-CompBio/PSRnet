setwd('~/project/PSVnet/V8/regulon/')
rm(list=ls())

#TF <- as.matrix(read.delim("tf_list.txt",header=F,sep="\t")[,2])
tf_file <- read.csv('~/project/PSVnet/tf_list.txt',sep = '\t',header = T)
load('1_data/3_all_net.RData')
load('1_data/5_target_based_data.RData')


##############  KB  ##############
# 
# ## universe
# 
# universe <- as.vector(unique(functional_KB_net$gene_id))
# universe.number <- length(unique(universe))
# 
# 
# ## hits
# 
# T3SS <- read.xlsx("5_funtion_filter/1_T3SS_1448A.xlsx", colNames = F)
# T3SS <- unique(T3SS$X4)
# 
# total.hits <- intersect(T3SS,universe)
# total.hits.number <- length(total.hits)
# 
# 
# ## master regulator
# 
# tf = unique(functional_KB_net$tf_name)
# mra_results <- c()
# i = 10
# for (i in 1:length(tf)){
#   print(tf[i])
#   #temp <- tryCatch(read.table(paste0("functional_target_KB/",tf[i],".txt"),header=F),error=function(e) NULL)
#   temp <- functional_KB_net[functional_KB_net$tf_name == tf[i], c(1:4)]
#   temp <- distinct(temp)
#   
#   target.genes <- as.vector(unique(temp$gene_id))
#   target.genes.number <- length(target.genes)
#   observed.Hits <- length(intersect(target.genes,T3SS))
#   expected.Hits <- (total.hits.number/universe.number) * target.genes.number
#   pval <- phyper(observed.Hits - 1, m = total.hits.number, n = universe.number - total.hits.number, k = target.genes.number, lower.tail=F )
#   mra_results <- rbind(mra_results,data.frame(tf[i],universe.number,target.genes.number,total.hits.number,expected.Hits,observed.Hits,pval))
# }
# pvals.adj <- p.adjust(mra_results$pval, method="BH")
# mra_results <- cbind(mra_results,pvals.adj)
# mra_results <- mra_results[order(mra_results[,7]),]
# 
# sig = as.character(mra_results[mra_results$pvals.adj < 0.05,1])
# KB_sub <- functional_KB_net[functional_KB_net$tf_name == sig,]
# 
# write.csv(mra_results, file="8_master_regulator/T6SS_mra_results_KB.csv")
# 
# KB_target <- merge(KB_sub[,c(3:4)], functional_KB_target[,c(1:3)], by.x = 'gene_id', by.y = 'target_id')
# KB_target = distinct(KB_target)
# 


##############  MM  ##############


temp = intersect(T3SS,as.vector(unique(functional_MM_net$tf_id)))
temp = union(temp,total.hits)



## universe

universe <- union(as.vector(unique(functional_MM_net$gene_id)),
                  as.vector(unique(functional_MM_net$tf_id)))
universe.number <- length(unique(universe))


## hits

T3SS <- read.xlsx("5_funtion_filter/1_T3SS_1448A.xlsx", colNames = F)
T3SS <- unique(T3SS$X4)

total.hits <- intersect(T3SS,universe)
total.hits.number <- length(total.hits)


## master regulator

tf = unique(functional_MM_net$tf_name)
mra_results <- c()
i = 10
for (i in 1:length(tf)){
  print(tf[i])
  #temp <- tryCatch(read.table(paste0("functional_target_MM/",tf[i],".txt"),header=F),error=function(e) NULL)
  temp <- functional_MM_net[functional_MM_net$tf_name == tf[i], c(1:4)]
  temp <- distinct(temp)
  
  target.genes <- as.vector(unique(temp$gene_id))
  target.genes.number <- length(target.genes)
  observed.Hits <- length(intersect(target.genes,T3SS))
  expected.Hits <- (total.hits.number/universe.number) * target.genes.number
  pval <- phyper(observed.Hits - 1, m = total.hits.number, n = universe.number - total.hits.number, k = target.genes.number, lower.tail=F )
  mra_results <- rbind(mra_results,data.frame(tf[i],universe.number,target.genes.number,total.hits.number,expected.Hits,observed.Hits,pval))
}
pvals.adj <- p.adjust(mra_results$pval, method="BH")
mra_results <- cbind(mra_results,pvals.adj)
mra_results <- mra_results[order(mra_results[,7]),]
mra_results = distinct(mra_results)
write.csv(mra_results, file="8_master_regulator/T6SS_mra_results_MM.csv")


sig = as.character(mra_results[mra_results$pvals.adj < 0.05,1])
MM_sub <- functional_MM_net[functional_MM_net$tf_name %in% sig,]

MM_target <- merge(MM_sub[,c(3:4)], functional_MM_target[,c(1:3)], by.x = 'gene_id', by.y = 'target_id')
MM_target = distinct(MM_target)



## T3SS 

T3SS_df <- read.xlsx("5_funtion_filter/1_T3SS_1448A.xlsx", colNames = F)

save(T3SS_df,
     KB_sub,MM_sub,
     KB_target,MM_target,
     file = '8_master_regulator/MRA.RData')




