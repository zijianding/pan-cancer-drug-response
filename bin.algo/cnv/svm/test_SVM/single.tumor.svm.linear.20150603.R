###functions###
bootstrap_sample <- function(obs,class,BS,method){
  #this function only generates BS samples
  #in each sample, two classes are balanced; or supply the minority class
  #the sample size equals to the number of obs(+1/-1 due to odd-even)
  #return a list
  #NOTICE: the observations can NOT be numeric!!!
  if( is.numeric(obs) )
  {
    print("The observation can NOT be numeric!")
  }
  
  resp = unique(class)
  
  pos_obs = obs[class==resp[2]]
  neg_obs = obs[class==resp[1]]
  
  if(method=="balance")
  {
    sample.num = round((length(pos_obs)+length(neg_obs))/2)
    
    bs_mat.sample = matrix(ncol=BS,nrow=(2*sample.num))
    bs_mat.resp = matrix(ncol=BS,nrow=(2*sample.num))
    
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
  if(method=="supply")
  {
    #sample the minority class, make it equal to the majority
    mino_class = resp[1] #neg_obs class
    if(length(neg_obs)>length(pos_obs))
    {
      mino_class = resp[2]
    }
    
    sample.num = length(obs) - 2*length(obs[class==mino_class])
    bs_mat.sample = matrix( ncol=BS,nrow=(length(obs)+sample.num) )
    bs_mat.resp = matrix(ncol=BS,nrow=length(obs)+sample.num)
    for(bs in 1:BS)
    {
      sample.pats = sample(obs[class==mino_class],size=sample.num,replace=T)
      sample.resp = rep(mino_class,times=sample.num)
      bs_mat.sample[,bs] = c(obs,sample.pats)
      bs_mat.resp[,bs] = c(class,sample.resp)
    }
    bs_mat = list(bs_mat.sample,bs_mat.resp)
    return(bs_mat)
    
  }
  
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

####load data###
#cluster
# args <- commandArgs(trailingOnly=TRUE)
# cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
# cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")
# test_fold = as.numeric(as.character(args[3]))#the current test fold, ranging from 1 to 10
# core.cancer = as.character(args[4])

#desktop
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv")
cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")
core.cancer = "CESC"

#both cluster and desktop
test_fold = test_fold + 3


#libraries
library(e1071)
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
models = seq(0.01,0.5,by=0.01)
costs = c(0.001,0.01,1,5,10,50,100,500,1000,5000)
param.mat = matrix(0,nrow=rep(models,each=length(costs)),ncol=rep(costs,times=length(models)))

folds = 5 # 5-fold CV for best model selection
BS = 100 #for bootstrap
best_param.mat = matrix(0,nrow=2,ncol=2)
dimnames(best_param.mat) = list(c("pr","roc"),c("cost","feature"))
score.mat = matrix(0,nrow=nrow(param.mat),ncol=length(train.pats)) 

for(fold in 1:5)
{
  
}


for(ca in 1:length(core.cancer))
{
  core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer[ca],]
  #average classification score for every observation obtained in cross validation
  score.mat = matrix(0,nrow=K,ncol=nrow(core.info)) 
  colnames(score.mat) = core.info$patient
  
  
  #core.res = list()
  #start cross validation
  for(cv in 4:8)
  {
    #find the test sample, and classes of each observation
    #data
    test.pats = as.character(core.info$patient[as.character(core.info[,cv])=="test"])
    test.dat = cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))]
    test.dat = test.dat[,test.pats]
    test.info = core.info[as.character(core.info[,cv])=="test",]
    test.mat = as.matrix(t(test.dat))
    #response
    test.resp = as.character(test.info$response[match(colnames(test.dat),as.character(test.info$patient))])
    test.resp.lab = vector(mode="numeric",length=length(test.resp))
    test.resp.lab[which(test.resp=="insensitive")] = 1
    test.resp.lab[which(test.resp=="sensitive")] = 0
        
    #find the train sample, and classes of each observation
    train.pats = as.character(core.info$patient[as.character(core.info[,cv])=="train"])
    train.info = core.info[as.character(core.info[,cv])=="train",]
    train.dat = cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))]
    train.dat = train.dat[,train.pats]
    
    #bootstrap sample
    pos_obs = as.character(train.info$patient[train.info$response=="insensitive"])
    neg_obs = as.character(train.info$patient[train.info$response=="sensitive"])
    sample.num = ncol(train.dat)
    bs_mat.sample = matrix("NULL",ncol=BS,nrow=(2*sample.num))
    bs_mat.resp = matrix("NULL",ncol=BS,nrow=(2*sample.num))
    for(bs in 1:BS)
    {
      #sample positive observations
      sample_pos = sample(pos_obs,size=sample.num,replace=T)
      sample_pos_resp = as.character( train.info$response[match(sample_pos,as.character(train.info$patient))] )
      #sample negative observations
      sample_neg = sample(neg_obs,size=sample.num,replace=T)
      sample_neg_resp = as.character( train.info$response[match(sample_neg,as.character(train.info$patient))] )
      
      bs_mat.sample[,bs] = c(sample_pos,sample_neg)
      bs_mat.resp[,bs] = c(sample_pos_resp,sample_neg_resp)
    }
    
    #calculate average classification score for the test sample
    all.score = array(rep(0,length(test.pats)*K*BS),c(length(test.pats),K,BS))
    for(bs in 1:BS)
    {
      train.dat = as.matrix(t(cisplatin.dat[,bs_mat.sample[,bs]]))
      #tmp.score = vector(mode="numeric",length=K)
      for(j in 1:K)
      {
        #tmp.df = data.frame(x=train.dat,y=as.factor(bs_mat.resp[,bs]))
        #svm.res = svm(y~.,data=tmp.df,kernel="linear",cost=costs[j],scale=FALSE,probability=TRUE)
        svm.res = svm(x=train.dat,y=as.factor(bs_mat.resp[,bs]),kernel="linear",cost=costs[j],scale=FALSE,probability=TRUE)
        #test.res = predict(object=svm.res,newx=test.mat,y=as.factor(test.resp))
        #tmp.df = data.frame(x=test.mat,y=as.factor(test.resp))
        test.res = predict(svm.res,newdata=test.mat,decision.values = TRUE,probability=TRUE)
        all.score[,j,bs] = attributes(test.res)$decision.values
      }
      #all.score[,,bs] = test.res
    }
    test.ix = match(rownames(test.mat),colnames(score.mat))
    for(i in 1:length(test.pats))
    {
      curr.mat = all.score[i,,]
      score.mat[,test.ix] = rowMeans(curr.mat)
    }
  }
  
  
  #all.res[[core.cancer[c]]] = score.mat
  pat_resp = as.character(core.info$response)
  pat_resp.lab = vector(length=length(pat_resp),mode="numeric")
  pat_resp.lab[which(pat_resp=="insensitive")] = 1
  pat_resp.lab[which(pat_resp=="sensitive")] = 0
  
  #choose best model by precision-recall curve
  auc.best = 0
  auc.best_ix = 1
  for(i in 1:K)
  {
    prauc = pr.curve(score.mat[i,],weights.class0=pat_resp.lab,curve=T)
    if( prauc$auc.integral > auc.best )
    {
      auc.best = prauc$auc.integral
      auc.best_ix = i
    }
  }
  prauc.best = pr.curve(score.mat[auc.best_ix,],weights.class0=pat_resp.lab,curve=T)
  tmp_str = paste("Single cancer",core.cancer[ca],"by SVM(linear) with cost",costs[auc.best_ix],sep=" ")
  plot(prauc.best,main=tmp_str)
  best.costs_pr[1,ca] = costs[auc.best_ix]
  
  #choose best model by roc curve
  auc.best = 0
  auc.best_ix = 1
  for(i in 1:K)
  {
    roc = roc.curve(score.mat[i,],weights.class0=pat_resp.lab,curve=T)
    if( roc$auc > auc.best )
    {
      auc.best = roc$auc
      auc.best_ix = i
    }
  }
  roc.best = roc.curve(score.mat[auc.best_ix,],weights.class0=pat_resp.lab,curve=T)
  tmp_str = paste("Single cancer",core.cancer[ca],"by SVM(linear) with cost", costs[auc.best_ix],sep=" ")
  plot(roc.best,main=tmp_str)
  best.costs_roc[1,ca] = costs[auc.best_ix]
  
}

pdf("performance.svm.linear.single_cancer.20150513.pdf")


#test each single cancer model on all other cancers
# other.info = cisplatin.info[cisplatin.info$cancer %in% other.cancer,]
# other.pats = as.character(other.info$patient)
# other.dat = cisplatin.dat[,match(other.pats,colnames(cisplatin.dat))]
# other.dat = other.dat[,other.pats]
# other.resp = as.character(other.info$response)
# other.resp_lab = vector(mode="numeric",length=length(other.resp))
# other.resp_lab[other.resp=="insensitive"] = 1
# other.resp_lab[other.resp=="sensitive"] = 0
# 
# for(c in 1:length(core.cancer))
# {
#   score.mat = matrix(0,nrow=BS,ncol=nrow(other.info))
#   
#   core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer[c],]
#   core.pats = as.character(core.info$patient)
#   core.resp = as.character(core.info$response)
#   core.dat = cisplatin.dat[,match(core.pats,colnames(cisplatin.dat))]
#   
#   
#   pos_obs = as.character(core.info$patient[core.info$response=="insensitive"])
#   neg_obs = as.character(core.info$patient[core.info$response=="sensitive"])
#   sample.num = ncol(core.dat)
#   bs_mat.sample = matrix("NULL",ncol=BS,nrow=(2*sample.num))
#   bs_mat.resp = matrix("NULL",ncol=BS,nrow=(2*sample.num))
#   for(bs in 1:BS)
#   {
#     #sample positive observations
#     sample_pos = sample(pos_obs,size=sample.num,replace=T)
#     sample_pos_resp = as.character( core.info$response[match(sample_pos,as.character(core.info$patient))] )
#     #sample negative observations
#     sample_neg = sample(neg_obs,size=sample.num,replace=T)
#     sample_neg_resp = as.character( core.info$response[match(sample_neg,as.character(core.info$patient))] )
#     
#     bs_mat.sample[,bs] = c(sample_pos,sample_neg)
#     bs_mat.resp[,bs] = c(sample_pos_resp,sample_neg_resp)
#   }
#   
#   for(bs in 1:BS)
#   {
#     train.dat = as.matrix(t(cisplatin.dat[,bs_mat.sample[,bs]]))
#     glm.res = glmnet(x=train.dat,y=as.factor(bs_mat.resp[,bs]),family="binomial",alpha=1,lambda=best.lambdas[c])
#     test.res = predict(object=glm.res,newx=t(other.dat),y=as.factor(other.resp),type="response")
#     score.mat[bs,] = test.res
#   }
#   
#   #draw precision-recall curve
#   score.final = colMeans(score.mat)
#   
#   prauc.best = pr.curve(score.final,weights.class0=other.resp_lab,curve=T)
#   tmp_str = paste("Test",core.cancer[c],"on other cancers with lasso logistic regression model",sep=" ")
#   plot(prauc.best,main=tmp_str)
# 
# }



dev.off()

#　From github 
