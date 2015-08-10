#setwd("C:/Users/zding/workspace/drug_sensitivity/data/omics.drug_centric/cnv")
args <- commandArgs(trailingOnly=TRUE)
###input data and information####
#cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
#cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")


core.cancer = c("LUAD","BLCA","MESO","STAD")
other.cancer = setdiff(as.character(unique(cisplatin.info$cancer)),core.cancer)
library(tree)
library(PRROC)

#how to store the data and plot graphs

#train tumor->test tumor -> train model parameters -> test tumor performance
#as list

#costs = c(0.001,0.01,1,5,10,50,100,500,1000,5000)
K = 1
#BS = 100

#all.res = list()
pdf("performance.tree.pancancer.pooled.baseline.20150514.pdf")

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
  #glm.res = glmnet(x=train.mat,y=as.factor(train.resp),family="binomial",alpha=1,lambda=lambdas)
  curr.test_res = matrix(0,ncol=length(test.pats),nrow=K)
  for(j in 1:K)
  {
    df.train = data.frame(train.mat,train.resp)
    fit.tree = tree(train.resp~.,df.train)
    #fit.cv_tree = cv.tree(fit.tree,FUN=prune.misclass)
    #min.error_ix = which.min(fit.cv_tree$dev)
    #min.best_size = fit.cv_tree$size[min.error_ix]
    #prune_fit.tree = prune.misclass(fit.tree,best=min.best_size)
    #test.res = predict(prune_fit.tree,data.frame(test.mat),type="vector")
    test.res = predict(fit.tree,data.frame(test.mat),type="vector")
    curr.test_res[j,] = test.res[,1]
  }

  test.ix = match(rownames(test.mat),colnames(score.mat))
  score.mat[,test.ix] = curr.test_res
  
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
tmp_str = "Pan-cancer(pooled data) by decision tree(PR curve best)"
plot(prauc.best,main=tmp_str)
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
tmp_str = "Pan-cancer(pooled data) by decision tree(ROC curve best)"
plot(roc.best,main=tmp_str)
auc_roc.best_ix = auc.best_ix

#check the performance of the best model(AUC under PR curve) in each core cancer
for(ca in 1:length(core.cancer))
{
  curr.pats = as.character(cisplatin.info$patient[as.character(cisplatin.info$cancer)==core.cancer[ca]])
  curr.ix = match(curr.pats,colnames(score.mat))
  curr.score_mat =score.mat[,curr.ix]
  curr.resp_lab = pat_resp.lab[curr.ix]
  
  #precision-recall curve
  if( is.matrix(curr.score_mat) )
  {
    prauc = pr.curve(curr.score_mat[auc_pr.best_ix,],weights.class0=curr.resp_lab,curve=T)
    tmp_str = paste("Performance of best model on",core.cancer[ca],sep=" ")
    plot(prauc,main=tmp_str)
    
    roc = roc.curve(curr.score_mat[auc_roc.best_ix,],weights.class0=curr.resp_lab,curve=T)
    plot(roc,main=tmp_str)
  }
  if( !is.matrix(curr.score_mat) )
  {
    prauc = pr.curve(curr.score_mat,weights.class0=curr.resp_lab,curve=T)
    tmp_str = paste("Performance of best model on",core.cancer[ca],sep=" ")
    plot(prauc,main=tmp_str)
    
    roc = roc.curve(curr.score_mat,weights.class0=curr.resp_lab,curve=T)
    plot(roc,main=tmp_str)
  }
  
  #roc curve
  
  
}


#check the best performance in each core cancer
for(ca in 1:length(core.cancer))
{
  curr.pats = as.character(cisplatin.info$patient[as.character(cisplatin.info$cancer)==core.cancer[ca]])
  curr.ix = match(curr.pats,colnames(score.mat))
  curr.score_mat =score.mat[,curr.ix]
  curr.resp_lab = pat_resp.lab[curr.ix]
  
  #by precision-recall curve
  curr_auc.best = 0
  curr_auc.best_ix = 1
  if( is.matrix(curr.score_mat) )
  {
    for(i in 1:K)
    {
      prauc = pr.curve(curr.score_mat[i,],weights.class0=curr.resp_lab,curve=T)
      if( prauc$auc.integral > auc.best )
      {
        auc.best = prauc$auc.integral
        auc.best_ix = i
      }
    }
    prauc.best = pr.curve(curr.score_mat[curr_auc.best_ix,],weights.class0=curr.resp_lab,curve=T)
    tmp_str = "Test on core cancer by decision tree"
    plot(prauc.best,main=tmp_str)
  }

  if( !is.matrix(curr.score_mat) )
  {
    prauc.best = pr.curve(curr.score_mat,weights.class0=curr.resp_lab,curve=T)
    tmp_str = "Test on core cancer by decision tree"
    plot(prauc.best,main=tmp_str)
  }
  
  
  #by roc curve
  curr_auc.best = 0
  curr_auc.best_ix = 1
  if( is.matrix(curr.score_mat) )
  {
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
    tmp_str = "Test on core cancer by decision tree"
    plot(roc.best,main=tmp_str)
  }
  if( !is.matrix(curr.score_mat) )
  {
    roc.best = roc.curve(curr.score_mat,weights.class0=curr.resp_lab,curve=T)
    tmp_str = "Test on core cancer by decision tree"
    plot(roc.best,main=tmp_str)
  }
  
}


#check the performance of the best model on other cancer
#fit the best model
# core.pats = as.character(core.info$patient)
# core.dat = cisplatin.dat[,match(core.pats,colnames(cisplatin.dat))]
# core.resp = as.character(core.info$response)
# core.dat_mat = as.matrix(t(core.dat))
# best.model = glmnet(x=core.dat_mat,y=as.factor(core.resp),family="binomial",alpha=1,lambda=best.lambda)
# best.model_roc = glmnet(x=core.dat_mat,y=as.factor(core.resp),family="binomial",alpha=1,lambda=best.lambda_roc)
# 
# #test on all other cancer
# other.info = cisplatin.info[cisplatin.info$cancer %in% other.cancer,]
# other.pats = as.character(other.info$patient)
# other.dat = cisplatin.dat[,match(other.pats,colnames(cisplatin.dat))]
# other.dat = other.dat[,other.pats]
# other.dat_mat = as.matrix(t(other.dat))
# other.resp = as.character(other.info$response)
# other.resp_lab = vector(mode="numeric",length=length(other.resp))
# other.resp_lab[other.resp=="insensitive"] = 1
# other.resp_lab[other.resp=="sensitive"] = 0
# 
# test.res = predict(object=best.model,newx=other.dat_mat,y=as.factor(other.resp),type="response")
# test.res_roc = predict(object=best.model_roc,newx=other.dat_mat,y=as.factor(other.resp),type="response")
# 
# prauc = pr.curve(test.res,weights.class0=other.resp_lab,curve=T)
# plot(prauc,main="Performance of best model(PR curve) on other cancer")
# 
# roc = roc.curve(test.res_roc,weights.class0=other.resp_lab,curve=T)
# plot(roc,main="Performance of best model(ROC curve) on other cancer")

dev.off()
