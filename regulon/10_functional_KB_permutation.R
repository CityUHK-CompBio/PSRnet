library(openxlsx)

rm(list=ls())
setwd('~/project/PSVnet/V8/regulon/12_permutation/')



load('~/project/PSVnet/gff.RData')
load('../1_data/5_target_based_data.RData')
load('../1_data/4_tf_based_data.RData')



dis <- functional_KB_target
dis_sub <- dis[which(dis$count > 1),]

tf_targets_num <- functional_KB_tf$count
names(tf_targets_num) = functional_KB_tf$tf_id

rand_times <- 100000



permutate <- function(row) {
  target <- as.character(row[1])   # for each target with > 1 tf
  print(target)
  
  num_tf <- as.numeric(row[4])
  tfs <- as.vector(unlist(strsplit(as.character(row[5])," ")))            # its tf
  
  ## creat a df to save the results of each tf
  rand_sum <- matrix(0,nrow=rand_times,ncol=nrow(gff_m))
  colnames(rand_sum) <- gff_m$locus_tag
  
  for(j in 1 : length(tfs)){
    tf_id <- tfs[j]      # for each tf
    num_target <- tf_targets_num[tf_id]         # n target of that tf
    #print(tf_id)
    #print(num_target)
    rand <- t(rmultinom(n=rand_times, size = num_target, prob = rep(1,nrow(gff_m))))     # random x times, every time get n genes from all gene sets
    rand[which(rand > 0)] <- 1               ##?
    rand_sum <- rand_sum + rand              # add results of each tf
  }
  
  ## for every target, was randomed (number_of_tf times) multiply (each tf x times)
  rand_target <- rand_sum[,target]       # in x random times, how many times the target was included
  times <- length(which(rand_target >= num_tf))     # if a target was randomed more than number_of_target times, count it
  pval <- times/rand_times
  c(t(row),times,pval)
}

#system.time({ temp <- apply(dis_sub, 1, permutate)})

library(parallel)
detectCores()
no_cores <- detectCores() /4
cl <- makeCluster(no_cores)
clusterExport(cl, list('rand_times','gff_m','tf_targets_num'))
system.time({
  result <- parApply(cl=cl, dis_sub, 1, permutate)
})
stopCluster(cl)

result = t(result)
result = as.data.frame(as.matrix(result))


## adjust p

#result <- read.delim("functional_KB_target.txt",header=F,sep="\t")
colnames(result) = c(colnames(dis_sub),'times','pval')

result$adj.pval <- p.adjust(result$pval,method="BH")

write.table(result,"functional_KB_target.txt",row.names=F,col.names=T,quote=F,sep="\t")
write.xlsx(result, 'functional_KB_target.xlsx', col.names = T,row.names = F)


result <- read.delim("functional_KB_target.txt",header=T,sep="\t")



