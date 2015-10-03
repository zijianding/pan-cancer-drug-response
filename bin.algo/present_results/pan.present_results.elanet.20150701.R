###load data###
##lap-top##
# input_folder = "C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/pool/molecular"
# info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt"
# curr_title = "Pan-cancer Analysis based on CNV Only" 
# pdf_file = "pan.cnv_only.elanet.pdf"
# file_pattern = "pan.elanet.mid_res.test_[0-9]*.20150701.txt"
#cluster#
args <- commandArgs(trailingOnly=TRUE)
info_file = args[1]
#input
input_folder = args[2]
#file_pattern = as.character(args[3])
#output
output_folder = args[3]
create_folder = args[4]
# curr_title = args[6]
# pdf_file = args[7]
# roc_file = args[8]

file_pattern = "pan.elanet.mid_res.test_[0-9]*.20150701.txt"
curr_title = "roc"
pdf_file = "roc.pdf"
roc_file = "roc.txt"

##both##
source("source_all.R")
cisplatin.info = read.table(info_file,header=T,quote="",sep="\t")
pan.files = list.files(path=input_folder,full.names=T,pattern=file_pattern)
pan_score = lapply(pan.files,read_table)
filenames = lapply(pan.files,find_name)
names(pan_score) = filenames
###plot curves###
dir.create(file.path(output_folder,create_folder),showWarnings = F)
setwd(file.path(output_folder,create_folder))
pdf(pdf_file)
##ROC curves for each data split##
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

#average ROC curve#
roc = ave_roc(roc_res,type="vertical",curr_title,measure="sd")
title(NULL,"Vertical average;std bar")
roc = ave_roc(roc_res,type="vertical",curr_title,measure="quantile")
title(NULL,"Vertical average;quantile bar")
roc = ave_roc(roc_res,type="threshold",curr_title,measure="sd")
title(NULL,"Threshold average;std bar")

write.table(roc[[1]],roc_file,row.names=F,col.names=F,sep="\t",quote=F )

test_fit = auc_random(roc_res,test="ttest")
test_fit = auc_random(roc_res,test="ztest")
test_fit = auc_random(roc_res,test="wilcox")


##plot error rate with cutoff=0.5##
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

##output each roc auc##
auc_res = lapply( roc_res,function(x)return(trapz(x[,2],x[,1])) )
write.table(auc_res,"all_auc.txt",quote=F,sep="\t",row.names=F,col.names=F)
