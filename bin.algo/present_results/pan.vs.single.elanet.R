###library###
library(doParallel)
library(foreach)
library(pracma)
no_cores = detectCores()
###functions###


###load data###
##path##
#laptop#
# info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/miRNAseq/cisplatin.miRNAseq_fold_cv.mat.txt"
# 
# pan_path = "C:/Users/zding/workspace/projects/drug_sensitivity/results/cisplatin/miRNA/elastic_net/wilcox.no_cancer/performance/whole_genome.0.05/"
# sin_path = "C:/Users/zding/workspace/projects/drug_sensitivity/results/cisplatin/miRNA/elastic_net/wilcox.no_cancer/performance/cesc.whole_genome.0.05/"
# 
# core.cancer = "CESC"
# pdf_file = "tmp.pdf"
# drug = "cisplatin"
# data_type = "miRNA"
#cluster#
args <- commandArgs(trailingOnly=TRUE)
info_file = args[1]

pan_path = args[2]
sin_path = args[3]

core.cancer = args[4]
pdf_file = args[5]
drug = "cisplatin"
data_type = args[6]
pan_pattern = "pan.elanet.mid_res.test_[0-9]*.20150701.txt"
single_pattern = pan_pattern
source("source_all.R")

##specific functions##
error_calc <- function(data, info,core.cancer)
{
  #data as a list of matrice
  #each matrix: row as samples, column as true class, prediciton,etc
  #info: all information of patients by the drug
  
  #debug
  #data = pan_cancer_res ; info = cisplatin.info
  test_num = length(data)
  
  res_list = list()
  all_cancer = c()
  
  for(K in 1:test_num)
  {
    patients = rownames(data[[K]])
    cancers = as.character( info$cancer[ match(patients,as.character(info$patient)) ] )
    #cancer_type = unique( cancers )
    cancer_type = core.cancer
    
    #sen_vec = vector( mode="numeric",length=length(cancer_type) )
    #insen_vec = vector( mode="numeric", length=length(cancer_type) )
    
    ix = which( cancers==cancer_type )
    curr_dat = data[[K]][ix,]
    
    #specifity
    sen = sum(curr_dat$indicator[curr_dat$true=="sensitive"]==T)
    sen_all = sum(curr_dat$true=="sensitive")
    sen_wrong = sen/sen_all
    
    #sensitivity
    insen = sum(curr_dat$indicator[curr_dat$true=="insensitive"]==T)
    insen_all = sum(curr_dat$true=="insensitive")
    insen_wrong = insen/insen_all
    
    sen_vec = sen_wrong
    insen_vec = insen_wrong
    
    
    
    
    curr_res = cbind(sen_vec,insen_vec)
    rownames(curr_res) = cancer_type
    colnames(curr_res) = c("specificity","sensitivity")
    
    res_list[[K]] = curr_res
    all_cancer  = union(all_cancer,cancer_type)
    
    
  }
  
  #pool all results
  all_cancer = unique(all_cancer)
  final_res = matrix(NA,nrow=2*length(all_cancer),ncol=test_num)
  rownames(final_res) = rep(all_cancer,each=2)
  resp_vec = rep(c("Specificity","Sensitivity"),times=length(all_cancer))
  
  for(K in 1:test_num)
  {
    curr_res = res_list[[K]]
    cancer_type = rownames(curr_res)
    
    for( i in 1:nrow(curr_res) )
    {
      cancer_ix = match(cancer_type[i],rownames(final_res))
      final_res[cancer_ix,K] = curr_res[i,1]
      final_res[cancer_ix+1,K] = curr_res[i,2]
    }
  }
  
  final = data.frame(rownames(final_res),resp_vec,final_res[,1])
  colnames(final) = c("cancer","Type","value")
  for(K in 2:test_num)
  {
    tmp = data.frame(rownames(final_res),resp_vec,final_res[,K])
    colnames(tmp) = c("cancer","Type","value")
    final = rbind( final,tmp )
  }
  
  #boxplot
  library(ggplot2)
  ggplot( aes(y=value,x=cancer),data=final ) + 
    geom_boxplot(aes(fill=Type)) + 
    geom_point(position=position_dodge(width=0.75), aes(fill=Type)) +
    ggtitle("Performance in Single Cancer") +
    xlab("cancer types") +
    theme_set(theme_gray(base_size=25))
  
  return(final)
  
  
}

##both##
pan.files = list.files(path=pan_path,full.names=T,pattern=pan_pattern)
pan_score = lapply(pan.files,read_table)
filenames = lapply(pan.files,find_name)
names(pan_score) = filenames

sin.files = list.files(path=sin_path,full.names=T,pattern=single_pattern)
sin_score = lapply(sin.files,read_table)
filenames = lapply(sin.files,find_name)
names(sin_score) = filenames

cisplatin.info = read.table(info_file,header=T,quote="",sep="\t")

core.pats = as.character(cisplatin.info$patient[cisplatin.info$cancer == core.cancer])

###plot curves###
pdf(pdf_file)
##pan-cancer on core cancer##
cl = makeCluster(no_cores)
registerDoParallel(cl)
pan_roc_res <- foreach( i=1:length(pan_score) ) %dopar%
{
  #for glmnet, second class is target class
  curr_pats_all = colnames(pan_score[[i]])
  curr_pats_core = intersect(curr_pats_all,core.pats)
  curr_pats_core.ix = match(curr_pats_core,curr_pats_all)
  
  pan_score[[i]] = pan_score[[i]][ , curr_pats_core.ix]
  curr_pats = colnames(pan_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  #roc_res[[i]] = ensemble_roc(pan_score[[i]],resp,"sensitive")
  return(ensemble_roc(pan_score[[i]],resp,"sensitive"))
}
stopImplicitCluster()
stopCluster(cl)
pan_on_single.roc = ave_roc(  pan_roc_res,type="vertical",measure="sd",
                              paste("Performance of Pan-cancer model on ",core.cancer,sep="") )

pan_auc = lapply(pan_roc_res,function(x){xval=x[,2];
                                         yval=x[,1];
                                         return(trapz(xval,yval))})
pan_auc = unlist(pan_auc)
pan_vec = rep("pan",length(pan_auc))
##single cancer##
cl = makeCluster(no_cores)
registerDoParallel(cl)
sin_roc_res <- foreach( i=1:length(sin_score) ) %dopar%
{
  #for glmnet, second class is target class
  curr_pats = colnames(sin_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  # roc_res[[i]] = ensemble_roc(sin_score[[i]],resp,"sensitive")
  return( ensemble_roc(sin_score[[i]],resp,"sensitive") )
}
stopImplicitCluster()
stopCluster(cl)
single.roc = ave_roc(sin_roc_res,type="vertical",measure="sd",
                       paste("Performance of single-cancer model on ",core.cancer,sep="") )
sin_auc = lapply(sin_roc_res,function(x){xval=x[,2]
                                         yval=x[,1]
                                         return(trapz(xval,yval))})
sin_auc = unlist(sin_auc)
sin_vec = rep("sin",length(sin_auc))
vec1 = c(pan_vec,sin_vec)
vec2 = c(pan_auc,sin_auc)
df_auc = data.frame(model=vec1,auc=vec2)
fig = ggplot(data=df_auc,aes(x=model,y=auc)) + 
  geom_boxplot(notch=T) + 
  labs(title=paste("Average ROC comparison of ",drug,"by",data_type,sep=" "),x="cancer data",y="AUC")+ 
  theme(text=element_text(size=15),
        plot.title=element_text(size=18,face="bold")) 
print(fig)
###plot curves###
#plot pan cancer roc vs single cancer# 
pan_xy = pan_on_single.roc[[1]]
sin_xy = single.roc[[1]]
plot(pan_xy[,2],pan_xy[,1],col="blue",
     xlim=c(0,1),ylim=c(0,1),xlab="FPR",ylab="TPR",
     main= NULL,"l" )
title(paste("Pan-cancer VS ",core.cancer,sep=""), "Vertial average")
lines(sin_xy[,2],sin_xy[,1],col="red","l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,col="gray")
legend("topleft",legend=c("pan-cancer",core.cancer,"random"),
       lty=c(1,1,2),col=c("blue","red","gray"))


##compare sensitivity and specific##
#pan-cancer results#
cl = makeCluster(no_cores)
registerDoParallel(cl)
pan_cancer_res<-foreach(i=1:length(pan_score)) %dopar%
{
  curr_pats = colnames(pan_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  #error_res[[i]] = identify_class(pan_score[[i]],resp,"sensitive")
  return( identify_class(pan_score[[i]],resp,"sensitive") )
}
stopImplicitCluster()
stopCluster(cl)
pan_error = error_calc(pan_cancer_res,cisplatin.info,core.cancer=core.cancer)

#single cancer resutls#
cl = makeCluster(no_cores)
registerDoParallel(cl)
single_cancer_res <- foreach(i=1:length(sin_score)) %dopar%
{
  curr_pats = colnames(sin_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  #error_res[[i]] = identify_class(sin_score[[i]],resp,"sensitive")
  return( identify_class(sin_score[[i]],resp,"sensitive") )
}
stopImplicitCluster()
stopCluster(cl)

sin_error = error_calc(single_cancer_res,cisplatin.info,core.cancer=core.cancer)

pan_error[,1] = "pan-cancer"
sin_error[,1] = "single-cancer"

all_error = rbind(pan_error,sin_error)
library(ggplot2)
ggplot( aes(y=value,x=cancer),data=all_error ) + 
  geom_boxplot(aes(fill=Type)) + 
  geom_point(position=position_dodge(width=0.75), aes(fill=Type)) +
  ggtitle(paste("Performance Comparison on ",core.cancer,sep="")) +
  xlab("cancer types") +
  theme_set(theme_gray(base_size=20))

dev.off()

