#setwd("C:/Users/zding/workspace/drug_sensitivity/data/omics.drug_centric/cnv")
args <- commandArgs(trailingOnly=TRUE)
###input data and information####
#cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
#cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")

#core.cancer = c("LUAD","BLCA","MESO","STAD")
core.cancer = c("CESC","LUAD","BLCA")
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
#BS = 100

#all.res = list()
pdf("prcurve.pancancer.pooled.baseline.20150510.pdf")

#core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer[c],]
core.info = cisplatin.info[as.character(cisplatin.info$cancer) %in% core.cancer,]
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
  train.resp = as.character(core.info$response[match(train.pats,as.character(core.info$patient))])
  
  train.mat = as.matrix(t(train.dat))
  glm.res = glmnet(x=train.mat,y=as.factor(train.resp),family="binomial",alpha=1,lambda=lambdas)
  test.res = predict(object=glm.res,newx=test.mat,y=as.factor(test.resp),type="response")
  
  test.ix = match(rownames(test.mat),colnames(score.mat))
  for(i in 1:K)
  {
    score.mat[i,test.ix] = test.res[,i]
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
tmp_str = paste("Pan-cancer(pooled data) by lasso logistic regression with lambda",lambdas[auc.best_ix],sep=" ")
plot(prauc.best,main=tmp_str)
best.lambda = lambdas[auc.best_ix]
auc_pr.best_ix = auc.best_ix

#choose best model by roc curve
auc.best = 0
auc.best_ix = 1
for(i in 1:K)
{
  roc = roc.curve(score.mat[i,],weights.class0=pat_resp.lab,curve=T)
  if( roc$auc > auc.best )
  {
    auc.best = auc.best
    auc.best_ix = i
  }
}
roc.best = roc.curve(score.mat[auc.best_ix,],weights.class0=pat_resp.lab,curve=T)
tmp_str = paste("Pan-cancer(pooled data) by lasso logistic regression",lambdas[auc.best_ix],sep=" ")
plot(roc.best,main=tmp_str)
best.lambda_roc = lambdas[auc.best_ix]
auc_roc.best_ix = auc.best_ix

#check the performance of the best model in each core cancer
for(c in 1:length(core.cancer))
{
  curr.pats = as.character(cisplatin.info$patient[as.character(cisplatin.info$cancer)==core.cancer[c]])
  curr.ix = match(curr.pats,colnames(score.mat))
  curr.score_mat =score.mat[,curr.ix]
  curr.resp_lab = pat_resp.lab[curr.ix]
  
  #precision-recall curve
  prauc = pr.curve(curr.score_mat[auc_pr.best_ix,],weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Performance of best model on",core.cancer[c],sep=" ")
  plot(prauc,main=tmp_str)
  
  #roc curve
  roc = roc.curve(curr.score_mat[auc_roc.best_ix,],weights.class0=curr.resp_lab,curve=T)
  plot(roc,main=tmp_str)
}


#check the best performance in each core cancer
for(c in 1:length(core.cancer))
{
  curr.pats = as.character(cisplatin.info$patient[as.character(cisplatin.info$cancer)==core.cancer[c]])
  curr.ix = match(curr.pats,colnames(score.mat))
  curr.score_mat =score.mat[,curr.ix]
  curr.resp_lab = pat_resp.lab[curr.ix]
  
  #by precision-recall curve
  curr_auc.best = 0
  curr_auc.best_ix = 1
  for(i in 1:K)
  {
    prauc = pr.curve(curr.score_mat[i,],weights.class0=curr.resp_lab, curve=T)
    if(prauc$auc.integral > curr_auc.best)
    {
      curr_auc.best = prauc$auc.integral
      curr_auc.best_ix = i
    }
  }
  
  prauc.best = pr.curve(curr.score_mat[curr_auc.best_ix,],weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Test on core cancer",core.cancer[c]," by lasso logistic regression(",lambdas[curr_auc.best_ix],")",sep=" ")
  plot(prauc.best,main=tmp_str)
  
  #by roc curve
  curr_auc.best = 0
  curr_auc.best_ix = 1
  for(i in 1:K)
  {
    roc = roc.curve(curr.score_mat[i,],weights.class0=curr.resp_lab, curve=T)
    if( roc$auc > curr_auc.best )
    {
      curr_auc.best = roc$auc
      curr_auc.best_ix = i
    }
  }
  roc.best = roc.curve(curr.score_mat[curr_auc.best_ix,],weights.class0=curr.resp_lab,curve=T)
  tmp_str = paste("Test on core cancer",core.cancer[c]," by lasso logistic regression(",lambdas[curr_auc.best_ix],")",sep=" ")
  plot(roc.best,main=tmp_str)
}


#check the performance of the best model on other cancer
#fit the best model
core.pats = as.character(core.info$patient)
core.dat = cisplatin.dat[,match(core.pats,colnames(cisplatin.dat))]
core.resp = as.character(core.info$response)
core.dat_mat = as.matrix(t(core.dat))
best.model = glmnet(x=core.dat_mat,y=as.factor(core.resp),family="binomial",alpha=1,lambda=best.lambda)
best.model_roc = glmnet(x=core.dat_mat,y=as.factor(core.resp),family="binomial",alpha=1,lambda=best.lambda_roc)

#test on all other cancer
other.info = cisplatin.info[cisplatin.info$cancer %in% other.cancer,]
other.pats = as.character(other.info$patient)
other.dat = cisplatin.dat[,match(other.pats,colnames(cisplatin.dat))]
other.dat = other.dat[,other.pats]
other.dat_mat = as.matrix(t(other.dat))
other.resp = as.character(other.info$response)
other.resp_lab = vector(mode="numeric",length=length(other.resp))
other.resp_lab[other.resp=="insensitive"] = 1
other.resp_lab[other.resp=="sensitive"] = 0

test.res = predict(object=best.model,newx=other.dat_mat,y=as.factor(other.resp),type="response")
test.res_roc = predict(object=best.model_roc,newx=other.dat_mat,y=as.factor(other.resp),type="response")

prauc = pr.curve(test.res,weights.class0=other.resp_lab,curve=T)
plot(prauc,main="Performance of best model(PR curve) on other cancer")

roc = roc.curve(test.res_roc,weights.class0=other.resp_lab,curve=T)
plot(roc,main="Performance of best model(ROC curve) on other cancer")

dev.off()
