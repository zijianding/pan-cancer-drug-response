#change data storage format
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/data/")
##the omics data
library(data.table)
ic50_dat = fread("en_input_w5.csv")
write.table(ic50_dat,"en_input_w5.txt",sep="\t",quote=F,row.names=F,col.names=F)
ic50_dat = fread("en_input_w5.txt",header=T)
ic50_df = data.frame(ic50_dat)
write.table(ic50_df,"en_input_w5.txt",sep="\t",quote=F,row.names=F,col.names=T)
##the sensitivity data
cisplatin.dat = fread("sensitivity_data_for_drug_1005.csv",header=T)
View(cisplatin.dat)
cisplatin = data.frame(cancer=cisplatin.dat$"Cancer Type",
                       patient=cisplatin.dat$"Cell Line Name",
                       drug=cisplatin.dat$"Drug Name",
                       response=cisplatin.dat$"IC 50")
write.table(cisplatin,"cisplatin.ic50.txt",col.names=T,row.names=F,quote=F,sep="\t")

#generate train and test data

#include test with correlation test

#glmnet multillinear in bootstrap to identify features

#multiliear regression in ensemble model

