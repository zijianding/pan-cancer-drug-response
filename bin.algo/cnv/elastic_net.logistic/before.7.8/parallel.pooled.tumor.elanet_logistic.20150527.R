####for CNV data####
###system parameters###
options(digits = 15)

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
test_fold = as.numeric(as.character(args[3]))#the current test fold, ranging from 1 to 5
output_folder = args[4] #no "/" at the end
#desktop
# setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv")
# cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
# cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")
# test_fold=1
# output_folder = "."
#both
test_fold = test_fold + 3

#libraries
library(glmnet)
library(PRROC)
library(doParallel)
library(foreach)
no_cores = detectCores()

###preprocess data###
#core.info = cisplatin.info[as.character(cisplatin.info$cancer) %in% core.cancer,]
core.info = cisplatin.info
core.cancer = c("CESC","LUAD", "BLCA")

#add the cancer type to cisplatin.dat(rbind)
cancer_type = cisplatin.info$cancer[match(colnames(cisplatin.dat),as.character(cisplatin.info$patient))]
tmp_df = data.frame(cancer_type,rep(1,times=length(cancer_type)))
colnames(tmp_df) = c("cancer","tmp")
cancers = t(model.matrix(tmp~cancer,tmp_df)[,-1])
colnames(cancers) = colnames(cisplatin.dat)
cisplatin.dat = rbind(cisplatin.dat,cancers)

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
BS = 100
folds = 5
null_model = 0.1
least_feat = 5
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

#average classification score for every observation obtained in cross validation
score.mat = matrix(0,nrow=nrow(mat.params),ncol=length(train.pats) )
colnames(score.mat) = train.pats
df.mat = matrix(NA,nrow=length(alphas),ncol=length(lambdas))
dimnames(df.mat) = list(alphas,lambdas)
df_freq.mat = matrix(NA,nrow=length(alphas),ncol=length(lambdas))
dimnames(df_freq.mat) = list(alphas,lambdas)

#data partition for CV
cv.mat = partition_data(train.pats,5)

for(fold in 1:folds)
{
  #print messegae
  tmp_str = paste("Current iteration:",fold,"of 5-fold cross validation",sep=" ")
  print(tmp_str)
  
  #find current validation data
  curr.validation_pats = rownames(cv.mat)[cv.mat[,fold]=="validation"]
  curr.validation_resp = as.character(train.info$response[match(curr.validation_pats,as.character(train.info$patient))])
  curr.validation_dat = as.matrix(t( cisplatin.dat[,curr.validation_pats] ))
  
  #bootstrap samples
  curr.train_pats = rownames(cv.mat)[cv.mat[,fold]=="train"]
  curr.train_resp = as.character(train.info$response[match(curr.train_pats,as.character(train.info$patient))])
  list.bs_mat = bootstrap_sample(curr.train_pats,curr.train_resp,BS=100)
  bs_mat.pats = list.bs_mat[[1]]
  bs_mat.resp = list.bs_mat[[2]]
  
               
  
  #a parallel version
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  curr_res <- foreach(alpha=rep(1:length(alphas),each=BS),
                      bs=rep(1:BS,times=length(alphas)),
                      .combine='comb', .multicombine=TRUE,
                      .init=list(list(), list()) ) %dopar%
  {
    library(glmnet)
    #find train dat
    curr.train_dat = as.matrix( t( cisplatin.dat[,match(bs_mat.pats[,bs],colnames(cisplatin.dat))] ) )
    
    #train model
    glm.res = glmnet(x=curr.train_dat,y=as.factor(bs_mat.resp[,bs]),family="binomial",
                     alpha=alphas[alpha],lambda=lambdas)
    #test model
    test.res = predict(object=glm.res,newx=curr.validation_dat,
                       y=as.factor(curr.validation_resp),type="response")
    
    if( levels(as.factor(bs_mat.resp[,bs]))[2] == "insensitive" )
    {
      #return(t(test.res))
      test.res = t(test.res)
    }
    if( levels(as.factor(bs_mat.resp[,bs]))[2] != "insensitive" )
    {
      #return(t(1-test.res))
      test.res = 1 - t(test.res)
    }
    #all.score[,((alpha-1)*K+1):(alpha*K),bs] = test.res
    return(list(test.res,glm.res$df))
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  list_score = curr_res[[1]]
  list_df = curr_res[[2]]
  
  all.score = array( data = NA, dim = c(nrow(mat.params),length(curr.validation_pats),BS) )
  all.df = array( data=NA,dim=c(length(alphas),length(lambdas),BS) )  
  
  for(i in 1:length(list_score) )
  {
    alpha_ix = ceiling(i/BS)
    bs = i%%BS
    if(bs==0)
    {
      bs=BS
      #all.score[,((alpha_ix-1)*K+1):(alpha_ix*K),bs] = curr_res[[i]]
      if( sum( dim(all.score[ ((alpha_ix-1)*K+1):(alpha_ix*K),,bs])==dim(list_score[[i]]) )==2 )
      {
        all.score[ ((alpha_ix-1)*K+1):(alpha_ix*K),,bs] = list_score[[i]]
        all.df[alpha_ix,,bs] = list_df[[i]]
      }
      if( sum( dim(all.score[ ((alpha_ix-1)*K+1):(alpha_ix*K),,bs])==dim(list_score[[i]]) )!=2 )
      {
        print("Wrong dimensions for assignment of parallel results!\n")
      }
      
    }
    if(bs!=0)
    {
      #all.score[,((alpha_ix-1)*K+1):(alpha_ix*K),bs] = curr_res[[i]]
      if( sum( dim(all.score[ ((alpha_ix-1)*K+1):(alpha_ix*K),,bs])==dim(list_score[[i]]) )==2 )
      {
        all.score[ ((alpha_ix-1)*K+1):(alpha_ix*K),,bs] = list_score[[i]]
        all.df[alpha_ix,,bs] = list_df[[i]]
      }
      if( sum( dim(all.score[ ((alpha_ix-1)*K+1):(alpha_ix*K),,bs])==dim(list_score[[i]]) )!=2 )
      {
        print("Wrong dimensions for assignment of parallel results!\n")
      }
      
    }
  }

  

  
  test.ix = match(curr.validation_pats,colnames(score.mat))
  for(i in 1:length(test.ix))
  {
    #curr.mat = all.score[i,,]
    curr.mat = all.score[,i,]
    #still the next codes are wrong....
    score.mat[,test.ix[i]] = rowMeans(curr.mat)
  }
  
  for( i in 1:length(alphas) )
  {
    #mean value of df
    curr_mat = all.df[i,,]
    df.mat[i,] = as.vector(rowMeans(curr_mat))
    
    #frequency of 0 df
    curr_mat[curr_mat!=0] = 1
    df_freq.mat[i,] = rowSums(curr_mat)
  }
}


#output classificatin scores and degree of freedom
#actually for debug
tmp_str = paste(output_folder,"/","pooled.test_",test_fold-3,".elanet.score_mat.CV.txt",sep="")
write.table(score.mat,tmp_str,col.names=T,row.names=T,quote=F,sep="\t")
# score.mat = read.table("../../../results/omics_features/cnv/cisplatin/elastic_net/pooled.test_1.elanet.score_mat(CV).txt",
#                        header=T,row.names=1,sep="\t",quote="")
# score.mat = as.matrix(score.mat)
#mean df under each parameter combination
tmp_str = paste(output_folder,"/","pooled.test_",test_fold-3,".elanet.degreeFreedom.CV.txt",sep="")
write.table(df.mat,tmp_str,col.names=T,row.names=T,quote=F,sep="\t")
df_freq.mat = 1 - df_freq.mat/BS
tmp_str = paste(output_folder,"/","pooled.test_",test_fold-3,".elanet.degreeFreedom_Freq.CV.txt",sep="")
write.table(df_freq.mat,tmp_str,col.names=T,row.names=T,quote=F,sep="\t")


#annotate drug response of the train data
pat_resp = as.character( core.info$response[match(colnames(score.mat),as.character(core.info$patient))] )
pat_resp.lab = vector(length=length(pat_resp),mode="numeric")
pat_resp.lab[which(pat_resp=="insensitive")] = 1
pat_resp.lab[which(pat_resp=="sensitive")] = 0



###identify the best model###
#record the best models(PR+ROC) CV classification score
best.classify_score_cv = matrix( 0,ncol=ncol(score.mat),nrow=2+2*length(core.cancer) )
rownames(best.classify_score_cv) = c("pool_pr",core.cancer,"pool_roc",core.cancer)
colnames(best.classify_score_cv) = colnames(score.mat)
#record parameters:whole, each caner
#PR, NOTICE: output t(best.param_pr)
best.param_pr = matrix( 0,nrow=length(core.cancer)+1,ncol=2 )
rownames(best.param_pr) = c("pool",core.cancer)
colnames(best.param_pr) = c("alpha","lambda")
#ROC, NOTICE: output t(best.param_roc)
best.param_roc = matrix( 0,nrow=length(core.cancer)+1,ncol=2 )
rownames(best.param_roc) = c("pool",core.cancer)
colnames(best.param_roc) = c("alpha","lambda")
#output the best curves of best models
file_name = paste(output_folder,"/","performance.elanet_logistic.pooled.test_",test_fold-3,".20150527.pdf",sep="")
pdf(file_name)

##best model performance of CV##
#PR curve#
auc.best = 0
auc.best_ix = 0
cl = makeCluster(no_cores)
registerDoParallel(cl)
auc.score<-foreach(i=1:nrow(mat.params),.combine='c') %dopar%
{
  library(PRROC)
  row_ix = match(mat.params[i,1],alphas)
  col_ix = match(mat.params[i,2],lambdas)
  if( (df.mat[row_ix,col_ix] <= least_feat) || (df_freq.mat[row_ix,col_ix] >= null_model) )
  {
    return(0)
  }
  if( (df.mat[row_ix,col_ix] > least_feat ) && (df_freq.mat[row_ix,col_ix] < null_model) )
  {
    prauc = pr.curve(score.mat[i,],weights.class0=pat_resp.lab,curve=T)
    return(prauc$auc.integral)
  } 
}
stopImplicitCluster()
stopCluster(cl)
auc.best = max(auc.score)
auc.best_ix = which( auc.score==max(auc.score) )
if( length(auc.best_ix>1) )
{
  auc.best_ix = auc.best_ix[length(auc.best_ix)]
}
best.param_pr[1,] = mat.params[auc.best_ix,]
prauc.best = pr.curve(score.mat[auc.best_ix,],weights.class0=pat_resp.lab,curve=T)
best.classify_score_cv[1,] = score.mat[auc.best_ix,]
tmp_str = "Pan-cancer(pooled data) by elasitic net logistic(CV)"
plot(prauc.best,main=tmp_str);
legend( "topright",legend=paste("alpha=",best.param_pr[1,1],"\tlambda=",best.param_pr[1,2],sep="") )

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
  if( (df.mat[row_ix,col_ix] <= least_feat) || (df_freq.mat[row_ix,col_ix]>=null_model) )
  {
    return(0)
  }
  if( (df.mat[row_ix,col_ix] > least_feat) && (df_freq.mat[row_ix,col_ix]<null_model) )
  {
    roc = roc.curve(score.mat[i,],weights.class0=pat_resp.lab,curve=T)
    return(roc$auc)
  }
}
stopImplicitCluster()
stopCluster(cl)
auc.best = max(auc.score)
auc.best_ix = which( auc.score==max(auc.score) )
if( length(auc.best_ix>1) )
{
  auc.best_ix = auc.best_ix[length(auc.best_ix)]
}
best.param_roc[1,] = mat.params[auc.best_ix,]
roc.best = roc.curve(score.mat[auc.best_ix,],weights.class0=pat_resp.lab,curve=T)
#这下边的index也是有问题，感觉应该是1+length(core.cancer)+1才对
#best.classify_score_cv[1+length(core.cancer),] = score.mat[auc.best_ix,]
best.classify_score_cv[1+length(core.cancer)+1,] = score.mat[auc.best_ix,]
tmp_str = "Pan-cancer by elasitic net logistic(CV)"
plot(roc.best,main=tmp_str);
legend( "bottomright",legend=paste("alpha=",best.param_roc[1,1],"\tlambda=",best.param_roc[1,2],sep="") )

#check the best model(whole) in each cancer
for(ca in 1:length(core.cancer))
{
  #find patients of current cancer
  curr.pats = train.pats[core.info$cancer[match(train.pats,as.character(core.info$patient))] == core.cancer[ca]]
  curr.ix = match(curr.pats,colnames(score.mat))
  curr.score = as.vector(best.classify_score_cv[1,curr.ix])
  curr.resp_lab = pat_resp.lab[curr.ix]
  
  #PR curve
  prauc = pr.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Performance of best model(PR) on",core.cancer[ca],sep=" ")
  plot(prauc,main=tmp_str);
  legend( "topright",legend=paste("alpha=",best.param_pr[1,1],"\tlambda=",best.param_pr[1,2],sep="") )
  
  #roc curve
  #这里可能有问题
  #curr.score = as.vector(best.classify_score_cv[1+length(core.cancer),curr.ix])
  curr.score = as.vector(best.classify_score_cv[1+length(core.cancer)+1,curr.ix])
  roc = roc.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Performance of best model(ROC) on",core.cancer[ca],sep=" ")
  plot(roc,main=tmp_str);
  legend( "bottomright",legend=paste("alpha=",best.param_roc[1,1],";lambda=",best.param_roc[1,2],sep="") )
 
}


#check the best performance in each cancer
for(ca in 1:length(core.cancer))
{
  #patients of current cancer and their classification scores
  curr.pats = train.pats[core.info$cancer[match(train.pats,as.character(core.info$patient))] == core.cancer[ca]]
  curr.ix = match(curr.pats,colnames(score.mat))
  curr.score_mat =score.mat[,curr.ix]
  curr.resp_lab = pat_resp.lab[curr.ix]
  
  #by precision-recall curve
  curr_auc.best = 0
  curr_auc.best_ix = 1
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  auc.score <- foreach(i=1:nrow(mat.params),.combine='c') %dopar%
  {
    library(PRROC)
    row_ix = match(mat.params[i,1],alphas)
    col_ix = match(mat.params[i,2],lambdas)
    if( (df.mat[row_ix,col_ix] <= least_feat) || (df_freq.mat[row_ix,col_ix]>=null_model) )
    {
      return(0)
    }
    if( (df.mat[row_ix,col_ix] > least_feat) && (df_freq.mat[row_ix,col_ix]<null_model) )
    {
      prauc = pr.curve( curr.score_mat[i,],weights.class0=curr.resp_lab,curve=T )
      return( prauc$auc.integral )
    }
  }
  stopImplicitCluster()
  stopCluster(cl)
  curr_auc.best = max(auc.score)
  curr_auc.best_ix = which(auc.score==max(curr_auc.best))
  if( length(curr_auc.best_ix)>1 )
  {
    curr_auc.best_ix = curr_auc.best_ix[length(curr_auc.best_ix)]
  }
  #plot PR curve
  prauc.best = pr.curve(curr.score_mat[curr_auc.best_ix,],weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Pan-cancer chosen by",core.cancer[ca],"(CV)",sep=" ")
  plot(prauc.best,main=tmp_str);
  best.param_pr[ca+1,] = mat.params[curr_auc.best_ix,]
  legend("topright",legend=paste("alpha=",best.param_pr[ca+1,1],"\tlambda=",best.param_pr[ca+1,2],sep=""))
  #record classification score
  best.classify_score_cv[ca+1,] = score.mat[curr_auc.best_ix,]
  
  #by roc curve
  curr_auc.best = 0
  curr_auc.best_ix = 1
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  auc.score <- foreach(i=1:nrow(mat.params),.combine='c') %dopar%
  {
    library(PRROC)
    row_ix = match(mat.params[i,1],alphas)
    col_ix = match(mat.params[i,2],lambdas)
    if( (df.mat[row_ix,col_ix] <= least_feat) || (df_freq.mat[row_ix,col_ix]>=null_model) )
    {
      return(0)
    }
    if( (df.mat[row_ix,col_ix] > least_feat) && (df_freq.mat[row_ix,col_ix]<null_model) )
    {
	  roc = roc.curve(curr.score_mat[i,],weights.class0=curr.resp_lab,curve=T)
      return(roc$auc)
    }
  }
  stopImplicitCluster()
  stopCluster(cl)
  curr_auc.best = max(auc.score)
  curr_auc.best_ix = which(auc.score==max(curr_auc.best))
  if( length(curr_auc.best_ix)>1 )
  {
    curr_auc.best_ix = curr_auc.best_ix[length(curr_auc.best_ix)]
  }
  best.param_roc[ca+1,] = mat.params[curr_auc.best_ix,]
  best.classify_score_cv[1+ca+length(core.cancer)+1,] = score.mat[curr_auc.best_ix,]
  roc.best = roc.curve(curr.score_mat[curr_auc.best_ix,],weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Pan-cancer(pooled) chosen by",core.cancer[ca],"(CV)",sep=" ")
  plot(roc.best,main=tmp_str);
  legend("bottomright",legend=paste("alpha=",best.param_roc[ca+1,1],"\tlambda=",best.param_roc[ca+1,2],sep=""))
  
}

#####best models on test data#####
print("Test best models")
#on the test cancer, actually 8 models here
list.bs_mat = bootstrap_sample(train.pats,train.resp,BS=100)
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]
#feature selection matrix
mat.feature_hitting = array(rep(0,nrow(train.dat)*100*8),c(nrow(train.dat),100,8))
rownames(mat.feature_hitting) = rownames(train.dat)

#directly testing on the test sample
test_score = array( NA, c(length(test.pats),100,8) )
for(bs in 1:BS)
{
  curr.train_pats = bs_mat.pats[,bs]
  curr.train_resp = bs_mat.resp[,bs]
  
  curr.train_dat = as.matrix( t( cisplatin.dat[,match(curr.train_pats,colnames(cisplatin.dat))] ) )
  #PR
  for(i in 1:(1+length(core.cancer)))
  {
    glm.res = glmnet(x=curr.train_dat,y=as.factor(curr.train_resp),family="binomial",
                     alpha=best.param_pr[i,1],lambda=best.param_pr[i,2])
    mat.feature_hitting[,bs,i] = as.vector(glm.res$beta)
    test.res = predict(object=glm.res,newx=test.mat,y=as.factor(test.resp),type="response")
    #all.score[,((alpha-1)*K+1):alpha*K,bs] = test.res
    if( levels(as.factor(curr.train_resp))[2]=="insensitive")
    {
      test_score[,bs,i] = test.res
    }
    if( levels(as.factor(curr.train_resp))[2]!="insensitive" )
    {
      test_score[,bs,i] = 1 - test.res
    }
    
  }
  #ROC
  for(i in 1:(1+length(core.cancer)))
  {
    glm.res = glmnet(x=curr.train_dat,y=as.factor(bs_mat.resp[,bs]),family="binomial",
                     alpha=best.param_roc[i,1],lambda=best.param_roc[i,2])
    mat.feature_hitting[,bs,i+1+length(core.cancer)] = as.vector(glm.res$beta)
    test.res = predict(object=glm.res,newx=test.mat,y=as.factor(test.resp),type="response")
    #all.score[,((alpha-1)*K+1):alpha*K,bs] = test.res
    if( levels(as.factor(curr.train_resp))[2] == "insensitive" )
    {
      test_score[,bs,i+1+length(core.cancer)] = test.res
    }
    if( levels(as.factor(curr.train_resp))[2] != "insensitive" )
    {
      test_score[,bs,i+1+length(core.cancer)] = 1 - test.res
    }
  }
    
}
test_score.ave = matrix(0,nrow=length(test.pats),ncol=2+2*length(core.cancer))
rownames(test_score.ave) = test.pats
colnames(test_score.ave) = c("pooled.pr",core.cancer,"pooled.roc",core.cancer)

for(i in 1:(2+2*length(core.cancer)))
{
  test_score.ave[,i] = rowMeans(test_score[,,i])
}
#characteristics of best model
mat.feature = matrix( 0,nrow=nrow(mat.feature_hitting),ncol=(2+2*length(core.cancer)) )
rownames(mat.feature) = rownames(train.dat)
colnames(mat.feature) = c("pooled.pr",core.cancer,"pooled.roc",core.cancer)
for( i in 1:(2+2*length(core.cancer)) )
{
  tmp_mat = mat.feature_hitting[,,i]
  tmp_mat[tmp_mat!=0] = 1
  mat.feature[,i] = rowSums(tmp_mat)/BS
}

###output test performance of best models###
#draw curves of two whole model
#下边这些地方都有问题啊，test_score.ave row是病人样本啊
#prauc.best = pr.curve(test_score.ave[1,],weights.class0=test.resp.lab,curve=T)
prauc.best = pr.curve(test_score.ave[,1],weights.class0=test.resp.lab,curve=T)
tmp_str = "Pan-cancer(pooled data) by elasitic net logistic(TEST)"
plot(prauc.best,main=tmp_str);
legend( "topright",legend=paste("alpha=",best.param_pr[1,1],";lambda=",best.param_roc[1,2],sep="") )

#roc.best = roc.curve(test_score.ave[2+length(core.cancer),],weights.class0=test.resp.lab,curve=T)
roc.best = roc.curve(test_score.ave[,2+length(core.cancer)],weights.class0=test.resp.lab,curve=T)
tmp_str = "Pan-cancer(pooled data) by elasitic net logistic(TEST)"
plot(roc.best,main=tmp_str);
legend( "bottomright",legend=paste("alpha=",best.param_roc[1,1],";lambda=",best.param_roc[1,2],sep="") )

#draw curves of each best model chosen by different core cancers
for(ca in 1:length(core.cancer))
{
  #on test data 
  curr.pats = test.pats[core.info$cancer[match(test.pats,as.character(core.info$patient))] == core.cancer[ca]]
  curr.ix = match(curr.pats,test.pats)
  curr.score_pr = test_score.ave[curr.ix,1+ca]
  curr.score_roc = test_score.ave[curr.ix,2+length(core.cancer)+ca]
  curr.resp_lab = test.resp.lab[curr.ix]
  
  #plot PR curve
  prauc.best = pr.curve(curr.score_pr,weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Pan-cancer chosen by",core.cancer[ca],"(TEST)",sep=" ")
  plot(prauc.best,main=tmp_str);
  legend("topright",legend=paste("alpha=",best.param_pr[ca+1,1],";lambda=",best.param_pr[ca+1,2],sep=""))

  
  #plot ROC curve
  roc.best = roc.curve(curr.score_roc,weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Pan-cancer chosen by",core.cancer[ca],"(TEST)",sep=" ")
  plot(roc.best,main=tmp_str);
  legend("bottomright",legend=paste("alpha=",best.param_roc[ca+1,1],";lambda=",best.param_roc[ca+1,2],sep=""))
}


###output model characteristics and classification scores###
#CV performance
file_name = paste(output_folder,"/","pan_pooled.elanet.cv.classify_score.test_",test_fold-3,".20150526.txt",sep="")
write.table(best.classify_score_cv,file=file_name,col.names=T,row.names=T,sep="\t",quote=F)

#test performance
file_name = paste(output_folder,"/","pan_pooled.elanet.test.classify_score.test_",test_fold-3,".20150526.txt",sep="")
write.table(test_score.ave,file=file_name,col.names=T,row.names=T,quote=F,sep="\t")

#features
file_name = paste(output_folder,"/","pan_pooled.elanet.feature_stat.test_",test_fold-3,".20150526.txt",sep="")
write.table(mat.feature,file=file_name,row.names=T,col.names=T,quote=F,sep="\t")

#best models
file_name = paste(output_folder,"/","pan_pooled.elanet.best_models.test_",test_fold-3,".20150526.txt",sep="")
best.param_pr = t(best.param_pr)
best.param_roc = t(best.param_roc)
best.param = cbind( best.param_pr, best.param_roc )
colnames(best.param) = c("pr",core.cancer,"roc",core.cancer)
write.table( best.param, file= file_name,row.names=T,col.names=T,quote=F,sep="\t" )

dev.off()
