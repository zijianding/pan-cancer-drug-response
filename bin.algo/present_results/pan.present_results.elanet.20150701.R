###library###
library(doParallel)
library(foreach)
no_cores = detectCores()
###functions###




###load data###
##lap-top##
# pan_path = "C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/pool/molecular"
# info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt"
# curr_title = "Pan-cancer Analysis based on CNV Only" 
# pdf_file = "pan.cnv_only.elanet.pdf"
# file_pattern = "pan.elanet.mid_res.test_[0-9]*.20150701.txt"
#cluster#
args <- commandArgs(trailingOnly=TRUE)
pan_path = args[1]
info_file = args[2]
curr_title = args[3]
pdf_file = args[4]
file_pattern = as.character(args[5])
roc_file = args[6]
##both##
cisplatin.info = read.table(info_file,header=T,quote="",sep="\t")
pan.files = list.files(path=pan_path,full.names=T,pattern=file_pattern)
pan_score = lapply(pan.files,read_table)
filenames = lapply(pan.files,find_name)
names(pan_score) = filenames
###plot curves###
pdf(pdf_file)
##ROC curve##
cl = makeCluster(no_cores)
registerDoParallel(cl)
roc_res <- foreach(i=1:length(pan_score) ) %dopar%
{
  #for glmnet, second class is target class
  curr_pats = colnames(pan_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  #roc_res[[i]] = ensemble_roc(pan_score[[i]],resp,"sensitive")
  return( ensemble_roc(pan_score[[i]],resp,"sensitive") )
}
stopImplicitCluster()
stopCluster(cl)

roc = ave_roc(roc_res,type="vertical",curr_title,measure="sd")
title(NULL,"Vertical average;std bar")
#for random comparison
write.table(roc[[1]],roc_file,row.names=F,col.names=F,sep="\t",quote=F )
#
roc = ave_roc(roc_res,type="vertical",curr_title,measure="quantile")
title(NULL,"Vertical average;quantile bar")
roc = ave_roc(roc_res,type="threshold",curr_title,measure="sd")
title(NULL,"Threshold average;std bar")

test_fit = auc_random(roc_res,test="ttest")
test_fit = auc_random(roc_res,test="ztest")
test_fit = auc_random(roc_res,test="wilcox")


##plot error rate##
cl = makeCluster(no_cores)
registerDoParallel(cl)
error_res <-foreach(i=1:length(pan_score)) %dopar%
{
  curr_pats = colnames(pan_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  #error_res[[i]] = identify_class(pan_score[[i]],resp,"sensitive")
  return( identify_class(pan_score[[i]],resp,"sensitive") )
}
stopImplicitCluster()
stopCluster(cl)
error_calc(error_res,cisplatin.info)

dev.off()
