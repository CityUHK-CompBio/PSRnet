library(varhandle)
library(openxlsx)

rm(list = ls())
setwd("~/project/PSVnet/V7/regulon/5_funtion_filter/")


######################## melt GO KEGG term ######################## 

## GO
go <- read.delim("1_gene_ontology_tab.txt",header=T,sep="\t")
go_accession <- unique(unfactor(go$Accession))

stat = NULL
accession <- go_accession[1]
go_gene <- list()
for (accession in go_accession){
  gene <- go[go$Accession == accession, "Locus.Tag"]
  gene <- unique(unfactor(gene))
  count <- length(gene)
  gene <- paste(gene, collapse = ' ')
  
  term <- go[go$Accession == accession, "GO.Term"]
  term <- unique(unfactor(term))
  
  stat = rbind(stat, data.frame(accession,gene,count,term))
  go_gene[term] = gene
}

write.table(stat,file = '2_GO_PS1448A',quote = F,row.names = F,col.names = F,sep = '\t')


## KEGG
kegg <- read.delim("1_kegg_pathways.txt",header=T,sep="\t")
kegg_accession <- unique(unfactor(kegg$Pathway.Xref))

stat = NULL
accession <- kegg_accession[1]
kegg_gene <- list()
for (accession in kegg_accession){
  gene <- kegg[kegg$Pathway.Xref == accession, "Locus.Tag"]
  gene <- unique(unfactor(gene))
  count <- length(gene)
  gene <- paste(gene, collapse = ' ')
  
  term <- kegg[kegg$Pathway.Xref == accession, "Pathway.Name"]
  term <- unique(unfactor(term))
  
  stat = rbind(stat, data.frame(accession,gene,count,term))
  kegg_gene[term] = gene
}

write.table(stat,file = '2_KEGG_PS1448A',quote = F,row.names = F,col.names = F,sep = '\t')



## T3SS
T3SS <- read.xlsx(xlsxFile = '1_T3SS_1448A.xlsx',colNames = F)

accession <- 'T3SS'
gene <- unique(T3SS$X4)
count <- length(gene)
gene <- paste(gene, collapse = ' ')
term <- 'T3SS'
stat = data.frame(accession,gene,count,term)

write.table(stat,file = '2_T3SS_PS1448A',quote = F,row.names = F,col.names = F,sep = '\t')










######################## GO KEGG jaccard ######################## 

go_kegg_gene <- append(go_gene, kegg_gene)

go <- read.delim("2_GO_PS1448A",header=F,sep="\t")
go_term <- as.vector(go$V4)
kegg <- read.delim("2_KEGG_PS1448A",header=F,sep="\t")
kegg_term <- as.vector(kegg$V4)
term <- c(go_term,kegg_term)
combination <-as.data.frame(t(combn(term,2)))
dim(combination)

count_jaccard <- function(term) {
  gene_1 <- as.vector(unlist(strsplit(unlist(go_kegg_gene[term[1]]), " ")))
  gene_2 <- as.vector(unlist(strsplit(unlist(go_kegg_gene[term[2]]), " ")))
  # print(gene_1)
  # print(gene_2)
  # print(intersect(gene_1,gene_2))
  # print(union(gene_1,gene_2))
  jaccard <- length(intersect(gene_1,gene_2))/length(union(gene_1,gene_2))
  jaccard
}
#term <- c('sequence-specific DNA binding','DNA replication')
test <- as.data.frame(apply(combination, 1, count_jaccard))
all_jaccard <- cbind(combination,test)

write.table(all_jaccard,file = '3_all_jaccard.txt',quote = F,row.names = F,col.names = F,sep = '\t')


save(all_jaccard,go_kegg_gene,kegg_gene,go_gene, file = 'modify_GO_KEGG.RData')



