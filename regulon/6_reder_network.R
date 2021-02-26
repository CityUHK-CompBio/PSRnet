rm(list=ls())
setwd('~/project/PSVnet/V8/regulon/1_data/')

load('3_all_net.RData')
load('5_target_based_data.RData')
load('5_target_based_data.RData')
load('3_functional_net.RData')

T3SS <- read.xlsx(xlsxFile = '../5_funtion_filter/1_T3SS_1448A.xlsx',colNames = F)
T3SS <- T3SS$X8

save(known,known_KB,known_MM,
     functional_KB_net,functional_MM_net,
     functional_KB_target,functional_MM_target,
     T3SS,
     file = '6_reder_network.RData')
