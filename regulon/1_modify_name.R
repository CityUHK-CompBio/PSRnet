
rm(list=ls())
setwd('~/project/PSVnet/V8/regulon/1_data/')

######### chipseq

cp chip_target from chip to regulon
remove gltR files
rename 0857 to tsiR



######### rna seq

cd V7/regulon/1_data
cp -r ../../rna/DEG* .

## DEG_KB

cd DEG_KB
mv 0857_0858_DEGs_adjp_with_log2FC.txt tsiS_DEGs_adjp_with_log2FC.txt
mv cvsR_cvsS_DEGs_adjp_with_log2FC.txt cvsS_DEGs_adjp_with_log2FC.txt
mv gacA_gacA_DEGs_adjp_with_log2FC.txt gacA_DEGs_adjp_with_log2FC.txt
mv gacA_gacS_DEGs_adjp_with_log2FC.txt gacS_DEGs_adjp_with_log2FC.txt
mv gltR_gtrS_DEGs_adjp_with_log2FC.txt ~/trash/
mv hrpS_2_DEGs_adjp_with_log2FC.txt hrpS_DEGs_adjp_with_log2FC.txt
mv ompR_ompR_DEGs_adjp_with_log2FC.txt ompR_DEGs_adjp_with_log2FC.txt
mv ompR_envZ_DEGs_adjp_with_log2FC.txt envZ_DEGs_adjp_with_log2FC.txt
mv phoP_phoQ_DEGs_adjp_with_log2FC.txt phoQ_DEGs_adjp_with_log2FC.txt
mv pilR_pilR_DEGs_adjp_with_log2FC.txt pilR_DEGs_adjp_with_log2FC.txt
mv pilR_pilS_DEGs_adjp_with_log2FC.txt pilS_DEGs_adjp_with_log2FC.txt
mv rhpR_rhpS_DEGs_adjp_with_log2FC.txt rhpS_DEGs_adjp_with_log2FC.txt



## DEG_MM

cd DEG_MM
mv 0857_0858_DEGs_adjp_with_log2FC.txt tsiS_DEGs_adjp_with_log2FC.txt
mv cvsR_cvsS_DEGs_adjp_with_log2FC.txt cvsS_DEGs_adjp_with_log2FC.txt
mv gacA_gacA_DEGs_adjp_with_log2FC.txt gacA_DEGs_adjp_with_log2FC.txt
mv gacA_gacS_DEGs_adjp_with_log2FC.txt gacS_DEGs_adjp_with_log2FC.txt
mv gltR_gtrS_DEGs_adjp_with_log2FC.txt ~/trash/
mv hrpS_2_DEGs_adjp_with_log2FC.txt hrpS_DEGs_adjp_with_log2FC.txt
mv ompR_ompR_DEGs_adjp_with_log2FC.txt ompR_DEGs_adjp_with_log2FC.txt
mv ompR_envZ_DEGs_adjp_with_log2FC.txt envZ_DEGs_adjp_with_log2FC.txt
mv phoP_phoQ_DEGs_adjp_with_log2FC.txt phoQ_DEGs_adjp_with_log2FC.txt
mv pilR_pilR_DEGs_adjp_with_log2FC.txt pilR_DEGs_adjp_with_log2FC.txt
mv pilR_pilS_DEGs_adjp_with_log2FC.txt pilS_DEGs_adjp_with_log2FC.txt
mv rhpR_rhpS_DEGs_adjp_with_log2FC.txt rhpS_DEGs_adjp_with_log2FC.txt





