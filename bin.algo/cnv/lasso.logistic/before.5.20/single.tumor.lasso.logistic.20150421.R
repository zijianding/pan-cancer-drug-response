#########################5.5###################################
#setwd("C:/Users/zding/workspace/drug_sensitivity/data/omics.drug_centric/cnv")
args <- commandArgs(trailingOnly=TRUE)
###input data and information####
#cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
#cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")


#core.cancer = c("LUAD","BLCA","MESO","STAD")
core.cancer = c("CESC","LUAD", "BLCA")
other.cancer = setdiff(as.character(unique(cisplatin.info$cancer)),core.cancer)
library(glmnet)
library(PRROC)

#how to store the data and plot graphs

#train tumor->test tumor -> train model parameters -> test tumor performance
#as list

lambdas = vector(mode="numeric",length=1001)
epsilon = 0.001
K = 1001
lambdas[K] = 0.5
lambdas[1] = 0.5*epsilon
for(i in 2:(K-1))
{
  lambdas[i] = 10^0.003 * lambdas[i-1]
}
BS = 100
best.lambdas = matrix(0,nrow=1,ncol=length(core.cancer))
best.lambdas_roc = matrix(0,nrow=1,ncol=length(core.cancer))
#colnames(best.lambdas) = best.lambdas

#all.res = list()
pdf("performance.logistic_lasso.single_cancer.20150519.pdf")
for(c in 1:length(core.cancer))
{
  core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer[c],]
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
    sample.num = round((length(pos_obs)+length(neg_obs))/2)
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
      glm.res = glmnet(x=train.dat,y=as.factor(bs_mat.resp[,bs]),family="binomial",alpha=1,lambda=lambdas)
      test.res = predict(object=glm.res,newx=test.mat,y=as.factor(test.resp),type="response")
      all.score[,,bs] = test.res
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
  tmp_str = paste("Single cancer",core.cancer[c],"by lasso logistic",lambdas[auc.best_ix],sep=" ")
  plot(prauc.best,main=tmp_str)
  best.lambdas[c] = lambdas[auc.best_ix]
  
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
  tmp_str = paste("Single cancer",core.cancer[c],"by lasso logistic", lambdas[auc.best_ix],sep=" ")
  plot(roc.best,main=tmp_str)
  best.lambdas_roc[c] = lambdas[auc.best_ix]
  
}




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
