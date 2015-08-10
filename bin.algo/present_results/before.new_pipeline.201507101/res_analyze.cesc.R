options(digits=15)
pdf("results_analsis.CESC.cisplatin_cnv_elanet.2015.6.22.pdf")
###cancer==CESC###
###functions###
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

###output line charts of AUC###
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
title("Peformance of Elastic Net Logistic Regression","Drug Cisplatin; Omics CNV; Response two levels")
#add a legend
legend(xrange[1],yrange[2],cex=0.8,col=colors,pch=plotchar,lty=linetype,
       title="Measurement",legend=colnames(mat_auc))

###output PR curves###
#CV performance
colors = rainbow(length(best_cvs))
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
  if(i==1)
  {
    plot(pr,color=colors[i],auc.main=F,lwd=0.5,main=NA,legend=F)
  }
  if(i>1)
  {
    plot(pr,color=colors[i],auc.main=F,lwd=0.5,add=T,legend=F)
  } 
}
title("PR curve of Elastic Net on CESC","Cisplatin;CNV;Cross Validation")
legend("topright",col=colors,legend=1:length(best_cvs),lty=rep(1,times=length(best_cvs)))


colors = rainbow(length(best_tests))
for(i in 1:length(best_tests))
{
  #test performance
  score_mat = best_tests[[i]]
  curr_ix = match(rownames(score_mat),as.character(cisplatin.info$patient))
  curr_info = cisplatin.info[curr_ix,]
  
  resp = vector(mode="numeric",length=nrow(curr_info))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  
  pr = pr.curve(score_mat[,1],weights.class0=resp,curve=T)
  if(i==1)
  {
    plot(pr,color=colors[i],auc.main=F,lwd=0.5,main=NA)
  }
  if(i>1)
  {
    plot(pr,color=colors[i],auc.main=F,lwd=0.5,add=T)
  }
  
  
}
title("PR curve of Elastic Net on CESC","Cisplatin;CNV;Test")
legend("topright",col=colors,legend=1:length(best_tests),lty=rep(1,times=length(best_tests)))

###output ROC curve###
colors = rainbow(length(best_cvs))
for(i in 1:length(best_cvs))
{
  #CV performance
  score_mat = best_cvs[[i]]
  curr_ix = match(rownames(score_mat),as.character(cisplatin.info$patient))
  curr_info = cisplatin.info[curr_ix,]
  
  resp = vector(mode="numeric",length=nrow(curr_info))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  
  roc = roc.curve(score_mat[,2],weights.class0=resp,curve=T)
  if(i==1)
  {
    plot(roc,color=colors[i],auc.main=F,lwd=0.5,main=NA)
  }
  if(i>1)
  {
    plot(roc,color=colors[i],auc.main=F,lwd=0.5,add=T)
  }
  
  
}
title("ROC curve of Elastic Net on CESC","Cisplatin;CNV;Cross Validation")
legend("bottomright",col=colors,legend=1:length(best_cvs),lty=rep(1,times=length(best_cvs)))

colors = rainbow(length(best_tests))
for(i in 1:length(best_tests))
{
  #CV performance
  score_mat = best_tests[[i]]
  curr_ix = match(rownames(score_mat),as.character(cisplatin.info$patient))
  curr_info = cisplatin.info[curr_ix,]
  
  resp = vector(mode="numeric",length=nrow(curr_info))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  
  roc = roc.curve(score_mat[,2],weights.class0=resp,curve=T)
  if(i==1)
  {
    plot(roc,color=colors[i],auc.main=F,lwd=0.5,main=NA)
  }
  if(i>1)
  {
    plot(roc,color=colors[i],auc.main=F,lwd=0.5,add=T)
  }
  
  
}
title("ROC curve of Elastic Net on CESC","Cisplatin;CNV;Test")
legend("bottomright",col=colors,legend=1:length(best_tests),lty=rep(1,times=length(best_tests)))

###output data and model characteristics###
lambdas = vector(mode="numeric",length=1001)
epsilon = 0.001
K = 1001
lambdas[K] = 0.5
lambdas[1] = 0.5*epsilon
for(i in 2:(K-1))
{
  lambdas[i] = 10^0.003 * lambdas[i-1]
}
lambdas = rev(lambdas)
alphas = seq(0,1,by=0.01)
res_mat = matrix(0,nrow=6,ncol=length(best_models))
dimnames(res_mat) = list(c("train","test","pr.df","pr.df.freq","roc.df","roc.df.freq"),
                         seq(1,length(best_models),by=1))
for(i in 1:length(best_models))
{
  curr_model = best_models[[i]]
  curr_df = df_list[[i]]
  dimnames(curr_df) = list(alphas,lambdas)
  curr.df_freq = df_freq_list[[i]]
  dimnames(curr.df_freq) = list(alphas,lambdas)
  
  ##data characteristics
  tmp = unlist(strsplit(names(best_models)[i],split="\\_"))
  tmp = tmp[length(tmp)]
  tmp = unlist(strsplit(tmp,split="\\."))
  curr_test = as.numeric(tmp[1])
  test_fold = 3 + curr_test
  #train
  score_mat = best_cvs[[i]]
  curr_ix = match(rownames(score_mat),as.character(cisplatin.info$patient))
  curr_info = cisplatin.info[curr_ix,]
  resp = vector(mode="numeric",length=nrow(curr_info))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  res_mat[1,i] = sum(resp==1)/length(resp)
  #test
  score_mat = best_tests[[i]]
  curr_ix = match(rownames(score_mat),as.character(cisplatin.info$patient))
  curr_info = cisplatin.info[curr_ix,]
  resp = vector(mode="numeric",length=nrow(curr_info))
  resp[as.character(curr_info$response)=="insensitive"] = 1
  resp[as.character(curr_info$response)!="insensitive"] = 0
  res_mat[2,i] = sum(resp==1)/length(resp)
  
  ##model characteristics
  #pr model
  alpha = curr_model[1,1]
  lambda = curr_model[1,2]
  row_ix = match(alpha,alphas)
  col_ix = match(round(lambda,digits=15),round(lambdas,digits=15))
  res_mat[3,i] = curr_df[row_ix,col_ix]
  res_mat[4,i] = curr.df_freq[row_ix,col_ix]
  #roc model
  alpha = curr_model[2,1]
  lambda = curr_model[2,2]
  row_ix = alpha*100+1
  col_ix = match(round(lambda,digits=15),round(lambdas,digits=15))
  res_mat[5,i] = curr_df[row_ix,col_ix]
  res_mat[6,i] = curr.df_freq[row_ix,col_ix] 
}

dev.off()
write.table(res_mat,"results_analysis.CESC.cisplatin_cnv_elanet.20150622.txt",
            col.names=T,row.names=T,quote=F,sep="\t")



