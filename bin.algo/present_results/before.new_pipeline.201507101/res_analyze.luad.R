
###cancer==LUAD###
###load data###
#read CV performance of CV
files = list.files(path="C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/LUAD/",
                   pattern="single.LUAD.elanet.cv_score.test_*",full.names=T)
read_table <- function(x)
{
  return( read.table(x,sep="\t",header=T,row.names=1,quote="") )
}
best_cvs = lapply(files,read_table)
find_name <- function(x)
{
  arr = unlist(strsplit(x,split="\\/",perl=T))
  return(arr[length(arr)])
}
filenames = lapply(files,find_name)
names(best_cvs) = filenames
#read test performance of best model
files = list.files(path="C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/LUAD/",
                   pattern="single.LUAD.elanet.test_score.test_*",full.names=TRUE)
best_tests = lapply(files,read_table)
filenames = lapply(files,find_name)
names(best_tests) = filenames
#read patient labels
patient.info = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt",
                          header=T,sep="\t",quote="")
cisplatin.info = patient.info[as.character(patient.info$cancer)=="LUAD",]
###calculate AUC###
library(PRROC)
mat_auc = matrix(0,nrow=length(best_cvs),ncol=4)
dimnames(mat_auc) = list(seq(1,length(best_cvs),by=1),c("pr_cv","pr_test","roc_cv","roc_test"))
for(i in 1:length(best_cvs))
{
  #CV performance
  score_mat = best_cvs[[i]]
  curr_ix = match(rownames(score_mat),as.character(cisplatin.info$patient))
  curr_info = cisplatin.info[curr_ix,]
  
  resp = vector(mode="numeric",length=nrow(curr_info))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  
  pr = pr.curve(score_mat[,1],weights.class0=resp,curve=T)
  roc = roc.curve(score_mat[,2],weights.class0=resp,curve=T)
  mat_auc[i,1] = pr$auc.integral
  mat_auc[i,3] = roc$auc
  
  #test performance
  score_mat = best_tests[[i]]
  curr_ix = match(rownames(score_mat),as.character(cisplatin.info$patient))
  curr_info = cisplatin.info[curr_ix,]
  
  resp = vector(mode="numeric",length=nrow(curr_info))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  
  pr = pr.curve(score_mat[,1],weights.class0=resp,curve=T)
  roc = roc.curve(score_mat[,2],weights.class0=resp,curve=T)
  mat_auc[i,2] = pr$auc.integral
  mat_auc[i,4] = roc$auc
  
}

###output figures###
xrange = range(0:(length(best_cvs)+1))
yrange = range(-0.001,1.001)
line_num = ncol(mat_auc)
plot(xrange,yrange,type="n",xlab="test set (Numbered)",ylab="AUC")
colors = rainbow(line_num)
linetype = c(1:line_num)
plotchar = seq(18,18+line_num,1)
#add lines
for(i in 1:line_num)
{
  lines(1:length(best_cvs), mat_auc[,i],type="b",lwd=1.5, lty=linetype[i],col=colors[i], pch=plotchar[i])
}
#add title
title("Peformance of Elastic Net Logistic Regression","Tumor LUAD;Drug Cisplatin; Omics CNV; Response two levels")
#add a legend
legend("bottomleft",cex=0.8,col=colors,pch=plotchar,lty=linetype,legend=colnames(mat_auc))
