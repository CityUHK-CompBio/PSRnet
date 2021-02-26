rm(list=ls())
setwd("~/project/PSVnet/V8/rna/")

load('../../rna_201906/mRNA_counts_201906.RData')
load('../../rna_201908/mRNA_counts_201908.RData')

load('../../gff.RData')

############## mrna counts KB ##############

# rna_201906 = ['hrpL','hrpR','hrpS_2','rhpR_rhpS','phoP_phoQ','mgrA','gacA_gacS','algU','pilR_pilR','pilR_pilS','ompR_ompR']
# rna_201908 = ['rpoN','0857_0858','ompR_envZ','gacA_gacA','psrA','aefR','vfr','cvsR_cvsS','gltR_gtrS']

# find names in identify_DEG.R
mRNA_counts.KB <- cbind(mRNA_counts_201906[,c("hrpL_KB_1","hrpL_KB_2","hrpR_KB_1","hrpR_KB_2",
                                              "psph_hrpS_MM_1","psph_hrpS_MM_2", "1448A_rhpS_KB_1","1448A_rhpS_KB_2",
                                              "psph_3729_KB_1","psph_3729_KB_2","mgrA_KB_1_18","mgrA_KB_1_4",
                                              "psph_3719_KB_1","psph_3719_KB_2","algU_KB_1","algU_KB_2",
                                              "pilR_KB_1","pilR_KB_2","190225_X495_FCHYC2HCCXY_L1_pilSKB_1","190225_X495_FCHYC2HCCXY_L1_pilSKB_2",
                                              "ompR_KB_1","ompR_KB_2")],
                     mRNA_counts_201908[,c("ox_rpoN_KB_1","ox_rpoN_KB_2","S0858_KB_1","S0858_KB_2",
                                           "envZ_KB_1","envZ_KB_2","gacA_KB_1","gacA_KB_2",
                                           "psrA_KB_1","psrA_KB_2","aefR_KB_1","aefR_KB_2",
                                           "Vfr_KB_1","Vfr_KB_2","cvsS_KB_1","cvsS_KB_2")])

temp <- merge(data.frame(id = rownames(mRNA_counts.KB)), gff_m[,c(10:12)], by.x = 'id', by = 'gene_ID')
rownames(mRNA_counts.KB) <- temp$locus_tag


############## mrna counts MM ##############

# rna_201906 = ['hrpL','hrpR','hrpS_2','rhpR_rhpS','phoP_phoQ','mgrA','gacA_gacS','algU','gacA_gacA']
# rna_201908 = ['rpoN','0857_0858','pilR_pilR','pilR_pilS','ompR_ompR','ompR_envZ','psrA','aefR','vfr','cvsR_cvsS','gltR_gtrS']

# find names in identify_DEG.R
mRNA_counts.MM <- cbind(mRNA_counts_201906[,c("hrpL_MM_1","hrpL_MM_2","hrpR_MM_1","hrpR_MM_2",
                                              "psph_hrpS_MM_1","psph_hrpS_MM_2","1448A_rhpS_MM_1","1448A_rhpS_MM_2",
                                              "psph_3729_MM_1","psph_3729_MM_2","mgrA_MM_1","mgrA_MM_2",
                                              "psph_3719_MM_1","psph_3719_MM_2","algU_MM_1","algU_MM_2",
                                              "psph_gacA_MM_1","psph_gacA_MM_2")],
                        mRNA_counts_201908[,c("ox_rpoN_MM_1","ox_rpoN_MM_2","S0858_MM_1","S0858_MM_2",
                                              "pilR_MM_1","pilR_MM_2","pilS_MM_1","pilS_MM_2",
                                              "ompR_MM_1","ompR_MM_2","envZ_MM_1","envZ_MM_2",
                                              "psrA_MM_1","psrA_MM_2","aefR_MM_1","aefR_MM_2",
                                              "Vfr_MM_1","Vfr_MM_2","cvsS_MM_1","cvsS_MM_2")])


temp <- merge(data.frame(id = rownames(mRNA_counts.MM)), gff_m[,c(10:12)], by.x = 'id', by = 'gene_ID')
rownames(mRNA_counts.MM) <- temp$locus_tag


save(mRNA_counts.KB,mRNA_counts.MM, 
     file = 'mRNA_counts.RData')
