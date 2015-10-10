pdf("pan_single_analysis.CESC.CNV.elanet.20150623.pdf")

###libraries###
library(PRROC)
###functions###
#options(digits=15)
read_table <- function(x)
{
  return( read.table(x,sep="\t",header=T,row.names=1,quote="",numerals="no.loss") )
}
find_name <- function(x)
{
  arr = unlist(strsplit(x,split="\\/",perl=T))
  return(arr[length(arr)])
}

###load data###
##pan-cancer##
pan_path = "C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/pooled/"
#cv#
pan_cv.files = list.files(path=pan_path,pattern="pan\\(pooled\\).elanet.cv.classify_score.test_[0-9].20150526", full.names=T)
pan_cv = lapply(pan_cv.files,read_table)
filenames = lapply(pan_cv.files,find_name)
names(pan_cv) = filenames
#test#
pan_test.files = list.files(path=pan_path,pattern="pan\\(pooled\\).elanet.test.classify_score.test_[0-9].20150526",full.names=T)
pan_test = lapply(pan_test.files,read_table)
filenames = lapply(pan_test.files,find_name)
names(pan_test) = filenames

##single cancer##
##CESC here##
#cv#
single_path = "C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/CESC/"
sgl_cv.files = list.files(single_path,pattern="single.CESC.elanet.cv_score.test_[0-9].20150609",full.names=T)
sgl_cv = lapply(sgl_cv.files,read_table)
filenames = lapply(sgl_cv.files,find_name)
names(sgl_cv) = filenames
#test#
sgl_test.files = list.files(single_path,pattern="single.CESC.elanet.test_score.test_[0-9].20150609",full.names=T)
sgl_test = lapply(sgl_test.files,read_table)
filenames = lapply(sgl_test.files,find_name)
names(sgl_test) = filenames

##data partition##
file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt"
patient.info = read.table(file,header=T,quote="",sep="\t")

###calculate AUC###
cancer = "CESC"
##pan-cancer on CESC##
#cv#
pan_cv.auc = matrix(0,ncol=2,nrow=length(pan_cv))
colnames(pan_cv.auc) = c("pr","roc")

pr_list = list()
roc_list = list()
for(i in 1:length(pan_cv))
{
  score_mat = pan_cv[[i]]
  tmp = unlist(strsplit(pan_cv.files[[i]],split="\\."))
  tmp = tmp[length(tmp)-2]
  tmp = unlist(strsplit(tmp,split="\\_"))
  tmp = tmp[length(tmp)]
  #find patients
  curr_info = patient.info[as.character(patient.info$cancer)==cancer,]
  test_fold = as.numeric(tmp)+3
  pats_ix = as.character(curr_info[,test_fold])=="train"
  curr_pats = as.character( curr_info$patient[pats_ix] )
  curr_info = curr_info[pats_ix,]
  row_ix = match(curr_pats,rownames(score_mat))
  score_mat = score_mat[row_ix,]
  
  resp = vector(mode="numeric",length=length(row_ix))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  
  #pr, col_ix = 2
  col_ix = 2
  pr = pr.curve(score_mat[,col_ix],weights.class0=resp,curve=T)
  pr_list[[i]] = pr
  pan_cv.auc[i,1] = pr$auc.integral

  #roc, col_ix = 6
  col_ix = 6
  roc = roc.curve(score_mat[,col_ix],weights.class0=resp,curve=T)
  roc_list[[i]] = roc
  pan_cv.auc[i,2] = roc$auc
}
#draw the curves
colors = rainbow(length(pr_list))
for(i in 1:length(pr_list))
{
  if(i==1)
  {
    plot(pr_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F)
  }
  if(i>1)
  {
    plot(pr_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F,add=T)
  }
}
title("PR curve of Elastic Net on CESC(Pan)","Cisplatin;CNV;Cross Validation")
legend("topright",col=colors,legend=1:length(pr_list),lty=rep(1,times=length(pr_list)))

colors = rainbow(length(roc_list))
for(i in 1:length(pr_list))
{
  if(i==1)
  {
    plot(roc_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F)
  }
  if(i>1)
  {
    plot(roc_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F,add=T)
  }
}
title("ROC curve of Elastic Net on CESC(Pan)","Cisplatin;CNV;Cross Validation")
legend("bottomright",col=colors,legend=1:length(pr_list),lty=rep(1,times=length(pr_list)))

#test#
pan_test.auc = matrix(0,ncol=2,nrow=length(pan_test))
colnames(pan_test.auc) = c("pr","roc")

pr_list = list()
roc_list = list()
for(i in 1:length(pan_test))
{
  score_mat = pan_test[[i]]
  tmp = unlist(strsplit(pan_test.files[[i]],split="\\."))
  tmp = tmp[length(tmp)-2]
  tmp = unlist(strsplit(tmp,split="\\_"))
  tmp = tmp[length(tmp)]
  #find patients
  curr_info = patient.info[as.character(patient.info$cancer)==cancer,]
  test_fold = as.numeric(tmp)+3
  pats_ix = as.character(curr_info[,test_fold])=="test"
  curr_pats = as.character( curr_info$patient[pats_ix] )
  curr_info = curr_info[pats_ix,]
  row_ix = match(curr_pats,rownames(score_mat))
  score_mat = score_mat[row_ix,]
  
  resp = vector(mode="numeric",length=length(row_ix))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  
  #pr, col_ix = 2
  col_ix = 2
  pr = pr.curve(score_mat[,col_ix],weights.class0=resp,curve=T)
  pr_list[[i]] = pr
  pan_test.auc[i,1] = pr$auc.integral
  
  #roc, col_ix = 6
  col_ix = 6
  roc = roc.curve(score_mat[,col_ix],weights.class0=resp,curve=T)
  roc_list[[i]] = roc
  pan_test.auc[i,2] = roc$auc
}
#draw the curves
colors = rainbow(length(pr_list))
for(i in 1:length(pr_list))
{
  if(i==1)
  {
    plot(pr_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F)
  }
  if(i>1)
  {
    plot(pr_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F,add=T)
  }
}
title("PR curve of Elastic Net on CESC(Pan)","Cisplatin;CNV;Test")
legend("topright",col=colors,legend=1:length(pr_list),lty=rep(1,times=length(pr_list)))

colors = rainbow(length(roc_list))
for(i in 1:length(pr_list))
{
  if(i==1)
  {
    plot(roc_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F)
  }
  if(i>1)
  {
    plot(roc_list[[i]],col=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F,add=T)
  }
}
title("ROC curve of Elastic Net on CESC(Pan)","Cisplatin;CNV;Test")
legend("bottomright",col=colors,legend=1:length(pr_list),lty=rep(1,times=length(pr_list)))

###load CESC data###
#read CV performance of CV
files = list.files(path="C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/CESC/",
                   pattern="single.CESC.elanet.cv_score.test_*",full.names=T)

best_cvs = lapply(files,read_table)
filenames = lapply(files,find_name)
names(best_cvs) = filenames
#read test performance of best model
files = list.files(path="C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/CESC/",
                   pattern="single.CESC.elanet.test_score.test_*",full.names=TRUE)
best_tests = lapply(files,read_table)
filenames = lapply(files,find_name)
names(best_tests) = filenames
#read patient labels
patient.info = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt",
                          header=T,sep="\t",quote="")
cisplatin.info = patient.info[as.character(patient.info$cancer)=="CESC",]
#read best parameters
files = list.files(path="C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/CESC/",
                   pattern="single.CESC.best_model.test_*",full.names=T)
best_models = lapply(files,read_table)
filenames = lapply(files,find_name)
names(best_models) = filenames
#read degree of freedom
files = list.files(path="C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/CESC/",
                   pattern="single.CESC.elanet.test_[0-9].elanet.degreeFreedom\\(CV\\)",full.names=T)
df_list = lapply(files,read_table)
filenames = lapply(files,find_name)
names(df_list) = filenames
#read frequency of degree of freedom
files = list.files(path="C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/elastic_net.logistic/single/CESC/",
                   pattern="single.CESC.elanet.test_[0-9].elanet.degreeFreedom_Freq",full.names=T)
df_freq_list = lapply(files,read_table)
filenames = lapply(files,find_name)
names(df_freq_list) = filenames

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


###Boxplot###
library(ggplot2)
cancer_type = vector(length=2*nrow(mat_auc)+nrow(pan_cv.auc)+nrow(pan_test.auc),
                     mode="character")
cancer_type[1:(2*nrow(mat_auc))] = "single"
cancer_type[-c(1:(2*nrow(mat_auc)))] = "pan"

test_type = vector(length=2*nrow(mat_auc)+nrow(pan_cv.auc)+nrow(pan_test.auc),
                   mode="character")
test_type[1:nrow(mat_auc)] = "CV"
test_type[(nrow(mat_auc)+1):(2*(nrow(mat_auc)))] = "Test"
test_type[(2*(nrow(mat_auc))+1):(2*(nrow(mat_auc))+nrow(pan_test.auc))] = "CV"
test_type[(2*(nrow(mat_auc))+nrow(pan_test.auc)+1):length(test_type)] = "Test"
val = vector(length=2*nrow(mat_auc)+nrow(pan_cv.auc)+nrow(pan_test.auc),
             mode="numeric")
val[1:nrow(mat_auc)] = mat_auc[,1]
val[(nrow(mat_auc)+1):(2*(nrow(mat_auc)))] = mat_auc[,2]
val[(2*(nrow(mat_auc))+1):(2*(nrow(mat_auc))+nrow(pan_test.auc))] = pan_cv.auc[,1]
val[(2*(nrow(mat_auc))+nrow(pan_test.auc)+1):length(test_type)] = pan_test.auc[,1]
pr.df = data.frame(cancer_type,test_type,val)
colnames(pr.df) = c("cancer","type","auc")
ggplot( aes(y=auc,x=type,fill=cancer),data=pr.df ) + 
  geom_boxplot() + 
  ggtitle("PR curve AUCs of Elastic Net logistic regression") +
  xlab("Cisplatin, CESC, CNV")






val = vector(length=2*nrow(mat_auc)+nrow(pan_cv.auc)+nrow(pan_test.auc),
             mode="numeric")
val[1:nrow(mat_auc)] = mat_auc[,3]
val[(nrow(mat_auc)+1):(2*(nrow(mat_auc)))] = mat_auc[,4]
val[(2*(nrow(mat_auc))+1):(2*(nrow(mat_auc))+nrow(pan_test.auc))] = pan_cv.auc[,2]
val[(2*(nrow(mat_auc))+nrow(pan_test.auc)+1):length(test_type)] = pan_test.auc[,2]
pr.df = data.frame(cancer_type,test_type,val)
colnames(pr.df) = c("cancer","type","auc")
ggplot( aes(y=auc,x=type,fill=cancer),data=pr.df ) + 
  geom_boxplot() + 
  ggtitle("ROC curve AUCs of Elastic Net logistic regression") +
  xlab("Cisplatin, CESC, CNV")

dev.off()