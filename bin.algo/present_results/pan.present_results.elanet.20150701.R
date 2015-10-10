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
#output
output_folder = args[3]
create_folder = args[4]
#pre-defined input and output files
file_pattern = "pan.elanet.mid_res.test_[0-9]*.20150701.txt"
param_file_pattern = "pan.elanet.param.test_[0-9]*.20150701.txt"
curr_title = "roc"
pdf_file = "roc.pdf"
roc_file = "roc.txt"
param_text = "param.txt"
source("source_all.R")

cisplatin.info = read.table(info_file,header=T,quote="",sep="\t")

#output open#
dir.create(file.path(output_folder,create_folder),showWarnings = F)
setwd(file.path(output_folder,create_folder))
pdf(pdf_file)

##auc distribution##
pan.files = list.files(path=input_folder,full.names=T,pattern=file_pattern)
pan_score = lapply(pan.files,read_table)
filenames = lapply(pan.files,find_name)
names(pan_score) = filenames

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
# roc = ave_roc(roc_res,type="vertical",curr_title,measure="quantile")
# title(NULL,"Vertical average;quantile bar")
# roc = ave_roc(roc_res,type="threshold",curr_title,measure="sd")
# title(NULL,"Threshold average;std bar")
write.table(roc[[1]],roc_file,row.names=F,col.names=F,sep="\t",quote=F )
# test_fit = auc_random(roc_res,test="ttest")
# test_fit = auc_random(roc_res,test="ztest")
# test_fit = auc_random(roc_res,test="wilcox")


##plot error rate with cutoff=0.5##
# cl = makeCluster(no_cores)
# registerDoParallel(cl)
# error_res <-foreach(i=1:length(pan_score)) %dopar%
# {
#   curr_pats = colnames(pan_score[[i]])
#   resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
#   #error_res[[i]] = identify_class(pan_score[[i]],resp,"sensitive")
#   return( identify_class(pan_score[[i]],resp,"sensitive") )
# }
# stopImplicitCluster()
# stopCluster(cl)
# error_calc(error_res,cisplatin.info)

dev.off()

##output each roc auc##
auc_res = lapply( roc_res,function(x)return(trapz(x[,2],x[,1])) )
write.table(auc_res,"all_auc.txt",quote=F,sep="\t",row.names=F,col.names=F)


#parameter collection#
dirs = list.dirs(path=input_folder,full.names=T,recursive=F)
output_names = lapply(dirs,
                      function(x){ tmp = strsplit(x,split="\\/"); 
                                   tmp=unlist(tmp); 
                                   return(tmp[length(tmp)]) } 
)
if(length(dirs)>0){
  for( j in 1:length(dirs))
  {
    curr_folder = dirs[j]
    output_file_name = output_names[[j]]
    
    pan.files = list.files(path=curr_folder,full.names=T,pattern=param_file_pattern)
    pan.files = sort(pan.files)
    param = lapply(pan.files,read_table)
    
    ##get paramters##
    p.thresh = unlist( lapply(param,function(x)return(x[1,2])) )
    sig.gene = unlist( lapply(param,function(x)return(x[2,2])) )
    thresh = unlist( lapply(param,function(x)return(x[3,2])) )
    freq.gene = unlist( lapply(param,function(x)return(x[4,2])) )
    
    
    df = data.frame(p.thresh=p.thresh,sig.gene=sig.gene,
                    thresh=thresh,freq.gene = freq.gene)
    
    write.table(df,output_file_name,row.names=F,col.names=T,sep="\t",quote=F)
  }
}else{
  pan.files = list.files(path=input_folder,full.names=T,pattern=param_file_pattern)
  pan.files = sort(pan.files)
  param = lapply(pan.files,read_table)
  
  ##get paramters##
  p.thresh = unlist( lapply(param,function(x)return(x[1,2])) )
  sig.gene = unlist( lapply(param,function(x)return(x[2,2])) )
  thresh = unlist( lapply(param,function(x)return(x[3,2])) )
  freq.gene = unlist( lapply(param,function(x)return(x[4,2])) )
  
  
  df = data.frame(p.thresh=p.thresh,sig.gene=sig.gene,
                  thresh=thresh,freq.gene = freq.gene)
  
  write.table(df,param_text,row.names=F,col.names=T,sep="\t",quote=F)
}

#combine the AUC and parameters#
df2 = data.frame(auc=unlist(auc_res),
                p.thresh=p.thresh,sig.gene=sig.gene,
                thresh=thresh,freq.gene = freq.gene)
write.table(df2,paste("all_res",nrow(df2),"txt",sep="."),row.names=F,col.names=T,sep="\t",quote=F)

