####functions###
bootstrap_sample <- function(obs,class,BS){
  #this function only generates BS samples
  #in each sample, two classes are balanced
  #the sample size equals to the number of obs(+1/-1 due to odd-even)
  #return a list
  resp = unique(class)
  
  pos_obs = obs[class==resp[2]]
  neg_obs = obs[class==resp[1]]
  sample.num = round((length(pos_obs)+length(neg_obs))/2)
  
  bs_mat.sample = matrix("NULL",ncol=BS,nrow=(2*sample.num))
  bs_mat.resp = matrix("NULL",ncol=BS,nrow=(2*sample.num))
  
  for(bs in 1:BS)
  {
    #sample positive observations
    sample_pos = sample(pos_obs,size=sample.num,replace=T)
    sample_pos_resp = class[match(sample_pos,obs)]
    #sample_pos_resp = as.character( train.info$response[match(sample_pos,as.character(train.info$patient))] )
    #sample negative observations
    sample_neg = sample(neg_obs,size=sample.num,replace=T)
    sample_neg_resp = class[match(sample_neg,obs)]
    #sample_neg_resp = as.character( train.info$response[match(sample_neg,as.character(train.info$patient))] )
    bs_mat.sample[,bs] = c(sample_pos,sample_neg)
    bs_mat.resp[,bs] = c(sample_pos_resp,sample_neg_resp)
  }
  
  bs_mat = list(bs_mat.sample,bs_mat.resp)
  return(bs_mat)
}

partition_data <- function(obs, k){
  #k-fold partition of obs
  cv.mat = matrix("NULL",nrow=length(obs),ncol=k)
  
  fd.size = floor(length(obs)/k)
  tmp = length(obs)%%k
  if(tmp>0)
  {
    fd.ix = rep(1:tmp,each = (fd.size+1))
    fd.ix = c(fd.ix,rep((tmp+1):k,each=fd.size))
  }
  if(tmp==0)
  {
    fd.ix = rep(1:k,each=fd.size)
  }
  
  for( j in 1:k)
  {
    cv.mat[fd.ix!=j,j] = "train"
    cv.mat[fd.ix==j,j] = "validation"
  }
  rownames(cv.mat) = obs
  return(cv.mat)
}

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

####load data####
#cluster
args <- commandArgs(trailingOnly=TRUE)
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")
test_fold = as.numeric(as.character(args[3]))
core.cancer = as.character(args[4])
output_folder = args[5] #no "/" at the end
#desktop
# setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv")
# cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
# cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")
# core.cancer = "BLCA"
# test_fold=1
#both
test_fold = test_fold + 3

#for debug
tmp_str = paste("The",test_fold,"test of Cancer",core.cancer,sep=" ")
print(tmp_str)
#

#libraries
library(glmnet)
library(PRROC)
library(doParallel)
library(foreach)
no_cores = detectCores()

###preprocess data###
other.cancer = setdiff(as.character(unique(cisplatin.info$cancer)),core.cancer)
core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer,]
#find the test sample, and classes of each observation
#test data
test.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="test"])
test.info = core.info[as.character(core.info[,test_fold])=="test",]
test.dat = cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))]
test.dat = test.dat[,test.pats]
test.mat = as.matrix(t(test.dat))
#test patients response
test.resp = as.character(test.info$response[match(colnames(test.dat),as.character(test.info$patient))])
test.resp.lab = vector(mode="numeric",length=length(test.resp))
test.resp.lab[which(test.resp=="insensitive")] = 1
test.resp.lab[which(test.resp=="sensitive")] = 0
#find the train sample, and classes of each observation
train.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="train"])
train.pats = sample(train.pats,size=length(train.pats),replace=FALSE)
train.info = core.info[as.character(core.info[,test_fold])=="train",]
train.dat = cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))]
train.dat = train.dat[,train.pats]
#responses
train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
train.resp.lab = vector(mode="numeric",length=length(train.resp))
train.resp.lab[which(train.resp=="insensitive")] = 1
train.resp.lab[which(train.resp=="sensitive")] = 0


###choose best model###
folds = 5 # 5-fold CV for best model selection
BS = 100 # 100 bootstrap samples one model
#elastic net models specified by parameters
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
mat.params = matrix(0,nrow=length(alphas)*K,ncol=2)
mat.params[,1] = rep(alphas,each=K)
mat.params[,2] = rep(lambdas,times=length(alphas))
#store the best model
best.param_pr = matrix(0,nrow=length(core.cancer),ncol=2) #only 1
best.param_roc = matrix(0,nrow=length(core.cancer),ncol=2) #only 1
#store the best characteristics of best model
file_name = paste(output_folder,"/","performance.elastic_net.logistic.single_",core.cancer,".test_",test_fold-3,".20150526.pdf",sep="")
pdf(file=file_name)
score.mat = matrix(0,nrow=nrow(mat.params),ncol=length(train.pats)) 
colnames(score.mat) = train.pats
#df.mat=average of degree of freedom under certain parameters
df.mat = matrix(NA,nrow=length(alphas),ncol=length(lambdas))
dimnames(df.mat) = list(alphas,lambdas)
#df_freq.mat= average frequency of 0 degree of freedom
df_freq.mat = matrix(NA,nrow=length(alphas),ncol=length(lambdas))
dimnames(df_freq.mat) = list(alphas,lambdas)
#data partition for CV
cv.mat = partition_data(train.pats,5)


for(fold in 1:folds)
{
  #print message
  tmp_str = paste("Current iteration:",fold,"of 5-fold cross validation",sep=" ")
  print(tmp_str)
  ###prepare data###
  #find current validation sets
  curr.validation_pats = rownames(cv.mat)[cv.mat[,fold]=="validation"]
  curr.validation_resp = as.character(train.info$response[match(curr.validation_pats,as.character(train.info$patient))])
  curr.validation_dat = as.matrix( t( cisplatin.dat[,match(curr.validation_pats,as.character(colnames(cisplatin.dat)))] ) )
  #find current train sets
  curr.train_pats = rownames(cv.mat)[cv.mat[,fold]=="train"]
  curr.train_resp = as.character(train.info$response[match(curr.train_pats,as.character(train.info$patient))])
  
  #bootstrap samples
  list.bs_mat = bootstrap_sample(curr.train_pats,curr.train_resp,BS=100)
  bs_mat.pats = list.bs_mat[[1]]
  bs_mat.resp = list.bs_mat[[2]]
  
  all.score = array( data=NA, dim=c(length(curr.validation_pats),nrow(mat.params),BS)  )
  all.df = array( data=NA,dim=c(length(alphas),length(lambdas),BS) )
  
  #here change to parallel version
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  curr_res <- foreach(alpha=rep(1:length(alphas),each=BS),
                      bs=rep(1:BS,times=length(alphas)),
                      .combine='comb', .multicombine=TRUE,
                      .init=list(list(), list()) ) %dopar%
  {
    library(glmnet)

    curr.train_dat = as.matrix( t( cisplatin.dat[,match(bs_mat.pats[,bs],colnames(cisplatin.dat))] ) )
    
    glm.res = glmnet(x=curr.train_dat,y=as.factor(bs_mat.resp[,bs]),family="binomial",
                     alpha=alphas[alpha],lambda=lambdas)
    
    test.res = predict(object=glm.res,newx=curr.validation_dat,
                       y=as.factor(curr.validation_resp),type="response")
    #all.score[,((alpha-1)*K+1):(alpha*K),bs] = test.res

    if( levels(as.factor(bs_mat.resp[,bs]))[2] != "insensitive"  )
    {
      test.res = 1 - test.res
    }
    return(list(test.res,glm.res$df))
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  list_score = curr_res[[1]]
  list_df = curr_res[[2]]
  
  for(i in 1:length(list_score) )
  {
    alpha_ix = ceiling(i/BS)
    bs = i%%BS
    if(bs==0)
    {
      bs=BS
      if( sum( dim(all.score[ ,((alpha_ix-1)*K+1):(alpha_ix*K),bs])==dim(list_score[[i]]) )==2 )
      {
        all.score[,((alpha_ix-1)*K+1):(alpha_ix*K),bs] = list_score[[i]]
        all.df[alpha_ix,,bs] = list_df[[i]]
      }
      if( sum( dim(all.score[ ,((alpha_ix-1)*K+1):(alpha_ix*K),bs])==dim(list_score[[i]]) )!=2 )
      {
        print("Wrong dimensions for assignment of parallel results!\n")
      }
    }
    if(bs!=0)
    {
      if( sum( dim(all.score[ ,((alpha_ix-1)*K+1):(alpha_ix*K),bs])==dim(list_score[[i]]) )==2 )
      {
        all.score[,((alpha_ix-1)*K+1):(alpha_ix*K),bs] = list_score[[i]]
        all.df[alpha_ix,,bs] = list_df[[i]]
      }
      if( sum( dim(all.score[ ,((alpha_ix-1)*K+1):(alpha_ix*K),bs])==dim(list_score[[i]]) )!=2 )
      {
        print("Wrong dimensions for assignment of parallel results!\n")
      }
    }
  }

  
  
  
  #average the output
  test.ix = match(curr.validation_pats,colnames(score.mat))
  for(i in 1:length(test.ix))
  {
    curr.mat = all.score[i,,]
    score.mat[,test.ix[i]] = rowMeans(curr.mat)
  }
  
  for( i in 1:length(alphas) )
  {
    #mean value of df
    curr_mat = all.df[i,,]
    df.mat[i,] = as.vector(rowMeans(curr_mat))
    
    #freq of NOT 0 for df
    curr_mat[curr_mat!=0] = 1
    df_freq.mat[i,] = rowSums(curr_mat)
  }
}

##CV classification score##
#actually for debug
#record the classification score of CV
tmp_str = paste(output_folder,"/","single.",core.cancer,".elanet.test_"
                ,test_fold-3,".score_mat(CV).20150609.txt",sep="")
write.table(score.mat,file=tmp_str,col.names=T,row.names=F,quote=F,sep="\t")
tmp_str = paste(output_folder,"/","single.",core.cancer,".elanet.test_"
                ,test_fold-3,".elanet.degreeFreedom(CV).txt",sep="")
write.table(df.mat,tmp_str,col.names=T,row.names=T,quote=F,sep="\t")
#Now df_freq.mat=freq of 0 for df
df_freq.mat = 1 - df_freq.mat/BS
tmp_str = paste(output_folder,"/","single.",core.cancer,".elanet.test_",
                test_fold-3,".elanet.degreeFreedom_Freq(CV).txt",sep="")
write.table(df_freq.mat,tmp_str,col.names=T,row.names=T,quote=F,sep="\t")


#annotate the train set label
pat_resp = as.character( core.info$response[match(train.pats,as.character(core.info$patient))] )
pat_resp.lab = vector(length=length(pat_resp),mode="numeric")
pat_resp.lab[which(pat_resp=="insensitive")] = 1
pat_resp.lab[which(pat_resp=="sensitive")] = 0

####best model performance of CV####
#PR curve
auc.best = 0
auc.best_ix = 0
#change to parallel
cl = makeCluster(no_cores)
registerDoParallel(cl)
auc.score<-foreach(i=1:nrow(mat.params),.combine='c') %dopar%
{
  library(PRROC)
  row_ix = match(mat.params[i,1],alphas)
  col_ix = match(mat.params[i,2],lambdas)
  if( (df.mat[row_ix,col_ix] <= 5) || (df_freq.mat[row_ix,col_ix]>=0.3) )
  {
    
    return(0)
  }
  if( (df.mat[row_ix,col_ix] > 5) && (df_freq.mat[row_ix,col_ix]<0.3) )
  {
    prauc = pr.curve(score.mat[i,],weights.class0=pat_resp.lab,curve=T)
    return(prauc$auc.integral)
  }
}
stopImplicitCluster()
stopCluster(cl)
auc.best = max(auc.score)
auc.best_ix = which.max(auc.score)
pr_best = as.vector(mat.params[auc.best_ix,])
prauc.best = pr.curve(score.mat[auc.best_ix,],weights.class0=pat_resp.lab,curve=T)
tmp_str = paste("Single cancer",core.cancer,"by elasitic net logistic(CV)",sep=" ")
plot(prauc.best,main=tmp_str);
legend( "topright",legend=paste("alpha=",pr_best[1],"\tlambda=",pr_best[2],sep="") )
#ROC curve
auc.best = 0
auc.best_ix = 1
cl = makeCluster(no_cores)
registerDoParallel(cl)
auc.score <-foreach(i=1:nrow(mat.params),.combine='c') %dopar%
{
  library(PRROC)
  row_ix = match(mat.params[i,1],alphas)
  col_ix = match(mat.params[i,2],lambdas)
  if( (df.mat[row_ix,col_ix] <= 1) || (df_freq.mat[row_ix,col_ix]>=0.7) )
  {
    return(0)
  }
  if( (df.mat[row_ix,col_ix] > 1) && (df_freq.mat[row_ix,col_ix]<0.7) )
  {
    roc = roc.curve(score.mat[i,],weights.class0=pat_resp.lab,curve=T)
    return(roc$auc)
  }
}
stopImplicitCluster()
stopCluster(cl)
auc.best = max(auc.score)
auc.best_ix = which.max(auc.score)
roc_best = as.vector(mat.params[auc.best_ix,])
roc.best = roc.curve(score.mat[auc.best_ix,],weights.class0=pat_resp.lab,curve=T)
tmp_str = paste("Single cancer",core.cancer,"by elasitic net logistic(CV)",sep=" ")
plot(roc.best,main=tmp_str);
legend( "bottomright",legend=paste("alpha=",roc_best[1],"\tlambda=",roc_best[2]) )



#####best models on test data#####
#best PR model
test_score.pr = matrix(0,nrow=length(test.pats),ncol=100)
mat_pr.feature_hitting = matrix(0,nrow=nrow(train.dat),ncol=100)
rownames(mat_pr.feature_hitting) = rownames(train.dat)
list.bs_mat = bootstrap_sample(train.pats,train.resp,BS=100)
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]

cl = makeCluster(no_cores)
registerDoParallel(cl)
test_res <- foreach( bs=1:100,.combine='comb', 
                    .multicombine=TRUE,
                    .init=list(list(), list()) ) %dopar%
{
  curr.train_pats = bs_mat.pats[,bs]
  curr.train_resp = bs_mat.resp[,bs]
  
  curr.train_dat = as.matrix( t( cisplatin.dat[,bs_mat.pats[,bs]] ) )
  library(glmnet)
  glm.res = glmnet(x=curr.train_dat,y=as.factor(bs_mat.resp[,bs]),family="binomial",
                   alpha=pr_best[1],lambda=pr_best[2])
  #mat_pr.feature_hitting[,bs] = as.vector(glm.res$beta)
  test.res = predict(object=glm.res,newx=test.mat,y=as.factor(test.resp),type="response")
  #all.score[,((alpha-1)*K+1):alpha*K,bs] = test.res
  #test_score.pr[,bs] = test.res
  if( levels(as.factor(bs_mat.pats[,bs]))[2] == "insensitive" )
  {
    
  }
  if( levels(as.factor(bs_mat.pats[,bs]))[2] != "insensitive" )
  {
    test.res = 1 - test.res
  }
  
  return(list(as.vector(glm.res$beta),test.res))
}
stopImplicitCluster()
stopCluster(cl)
list_pr.feature_hitting = test_res[[1]]
list_test_score = test_res[[2]]
for(bs in 1:100)
{
  mat_pr.feature_hitting[,bs] = list_pr.feature_hitting[[bs]]
  test_score.pr[,bs] = list_test_score[[bs]]
}
test_score.pr_ave = rowMeans(test_score.pr)
mat_pr.feature_hitting[mat_pr.feature_hitting!=0] = 1
mat_pr.feature = rowSums(mat_pr.feature_hitting)
mat_pr.feature = mat_pr.feature/BS
#plot PR curve
prauc.best = pr.curve(test_score.pr_ave,weights.class0=test.resp.lab,curve=T)
tmp_str = paste("Single cancer",core.cancer,"by elasitic net logistic(TEST)",sep=" ")
plot(prauc.best,main=tmp_str);
legend( "topright",legend=paste("alpha=",pr_best[1],"\tlambda=",pr_best[2],sep="") )


###best ROC model
test_score.roc = matrix(0,nrow=length(test.pats),ncol=100)
mat_roc.feature_hitting = matrix(0,nrow=nrow(train.dat),ncol=100)
rownames(mat_roc.feature_hitting) = rownames(train.dat)
list.bs_mat = bootstrap_sample(train.pats,train.resp,BS=100)
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]

cl = makeCluster(no_cores)
registerDoParallel(cl)
test_res <- foreach( bs=1:100,.combine='comb', 
                     .multicombine=TRUE,
                     .init=list(list(), list()) ) %dopar%
{
  curr.train_pats = bs_mat.pats[,bs]
  curr.train_resp = bs_mat.resp[,bs]
  
  curr.train_dat = as.matrix( t( cisplatin.dat[,bs_mat.pats[,bs]] ) )
  library(glmnet)
  ###########################################################################################
  glm.res = glmnet(x=curr.train_dat,y=as.factor(bs_mat.resp[,bs]),family="binomial",
                   alpha=roc_best[1],lambda=roc_best[2])
  ############################################################################################
  test.res = predict(object=glm.res,newx=test.mat,y=as.factor(test.resp),type="response")
  if( levels(as.factor(test.resp))[2] != "insensitive"  )
  {
    test.res = 1 - test.res
  }
  #mat_roc.feature_hitting[,bs] = as.vector(glm.res$beta)
  #all.score[,((alpha-1)*K+1):alpha*K,bs] = test.res
  #test_score.roc[,bs] = test.res
  return(list(as.vector(glm.res$beta),test.res))
}
stopImplicitCluster()
stopCluster(cl)
list_roc.feature_hitting = test_res[[1]]
list_test_score = test_res[[2]]
for(bs in 1:100)
{
  mat_roc.feature_hitting[,bs] = list_roc.feature_hitting[[bs]]
  test_score.roc[,bs] = list_test_score[[bs]]
}
test_score.roc_ave = rowMeans(test_score.roc)
mat_roc.feature_hitting[mat_roc.feature_hitting!=0] = 1
mat_roc.feature = rowSums(mat_roc.feature_hitting)
mat_roc.feature_hitting = mat_roc.feature_hitting/BS
#plot ROC curve
roc.best = roc.curve(test_score.roc_ave,weights.class0=test.resp.lab,curve=T)
tmp_str = paste("Single cancer",core.cancer,"by elasitic net logistic(TEST)",sep=" ")
plot(roc.best,main=tmp_str);
legend( "bottomright",legend=paste("alpha=",roc_best[1],"\tlambda=",roc_best[2]) )


dev.off()

###output model characteristics###
##CV best performance##
classify_score.best_pr = score.mat[auc.best_ix,]
classify_score.best_roc = score.mat[auc.best_ix,]
classify_score.best = rbind(classify_score.best_pr,classify_score.best_roc)
classify_score.best = t(classify_score.best)
colnames(classify_score.best) = c("pr","roc")
file_name = paste(output_folder,"/","single.",core.cancer,
                  ".elanet.cv_score.test_",test_fold-3,".20150609.txt",sep="")
write.table(classify_score.best,file=file_name,col.names=T,row.names=T,sep="\t",quote=F)

##test performance##
test_score.all = cbind(test_score.pr_ave,test_score.roc_ave)
colnames(test_score.all) = c("pr","roc")
rownames(test_score.all) = test.pats
file_name = paste(output_folder,"/","single.",core.cancer,
                  ".elanet.test_score.test_",test_fold-3,".20150609.txt",sep="")
write.table(test_score.all,file=file_name,col.names=T,row.names=T,quote=F,sep="\t")

#features
mat.feature = cbind(mat_pr.feature,mat_roc.feature)
colnames(mat.feature)=c("pr","roc")
file_name = paste(output_folder,"/","single.",core.cancer,
                  ".elanet.feature_stat.test_",test_fold-3,".20150609.txt",sep="")
write.table(mat.feature,file=file_name,row.names=T,col.names=T,quote=F,sep="\t")

#best models parameters
param_best.mat = rbind(pr_best,roc_best)
dimnames(param_best.mat) = list(c("pr","roc"),c("alpha","lambda"))
file_name = paste(output_folder,"/","single.",core.cancer,
                  ".elanet.best_parameters.test_",test_fold-3,".20150609.txt",sep="")



#　From github 
