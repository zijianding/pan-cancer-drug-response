setwd("C:/Users/zding/workspace/drug_sensitivity/data/omics.drug_centric/cnv")
###input data and information####
cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")

core.cancer = c("LUAD","BLCA","MESO","STAD")
library(glmnet)
library(PRROC)

#how to store the data and plot graphs

#train tumor->test tumor -> train model parameters -> test tumor performance
#as list

all.res = list()
for(c in 1:length(core.cancer))
{
  core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer[c],]
  core.res = list()
  for(cv in 4:8)
  {
    
    test.pats = as.character(core.info$patient[as.character(core.info[,cv])=="test"])
    test.dat = cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))]
    test.dat = test.dat[,test.pats]
    test.info = core.info[as.character(core.info[,cv])=="test",]
    
    test.resp = as.character(test.info$response[match(colnames(test.dat),as.character(test.info$patient))])
    test.resp.lab = vector(mode="numeric",length=length(test.resp))
    test.resp.lab[which(test.resp=="insensitive")] = 1
    test.resp.lab[which(test.resp=="sensitive")] = 0
    
    train.pats = as.character(core.info$patient[as.character(core.info[,cv])=="train"])
    train.info = core.info[as.character(core.info[,cv])=="train",]
    train.dat = cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))]
    train.dat = train.dat[,train.pats]
    
    #sampling the minority
    #pos.neg.ratio = sum(as.character(core.info$response)=="insensitive")/sum(as.character(core.info$response)=="sensitive")
    pos.neg.ratio = sum(as.character(train.info$response)=="insensitive")/sum(as.character(train.info$response)=="sensitive")
    if( pos.neg.ratio <= 0.6)
    {
      train.minority = as.character(train.info$patient[train.info$response=="insensitive"])
      #sample.num = nrow(train.info) - length(train.minority)
      sample.num = sum(as.character(train.info$response)=="sensitive") - sum(as.character(train.info$response)=="insensitive")
      sample.pats = sample(train.minority,size=sample.num,replace=T)
      col.num = match(sample.pats,colnames(train.dat))
      sample.dat = train.dat[,col.num]
      #train.dat = cbind(train.dat,sample.dat)
      train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
      sample.resp = as.character(train.info$response[match(sample.pats,as.character(train.info$patient))])
      train.resp = c(train.resp,sample.resp)
      train.dat = cbind(train.dat,sample.dat)
      #debug
      train.resp.lab = vector(mode="numeric",length=length(train.resp))
      train.resp.lab[which(train.resp=="insensitive")] = 1
      train.resp.lab[which(train.resp!="insensitive")] = 0
      tmp.char = paste("Current CV round",(cv-3),"of cancer",core.cancer[c],"with",sum(train.resp.lab),"positive samples",sep=" ")
      print(tmp.char)
      
    }
    if( pos.neg.ratio > 0.6)
    {
      train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
      #debug
      train.resp.lab = vector(mode="numeric",length=length(train.resp))
      train.resp.lab[which(train.resp=="insensitive")] = 1
      train.resp.lab[which(train.resp!="insensitive")] = 0
      tmp.char = paste("Current CV round",(cv-3),"of cancer",core.cancer[c],"with",sum(train.resp.lab),"positive samples",sep=" ")
      print(tmp.char)
    }
    
    
    
    train.mat = as.matrix(train.dat)
    test.mat = as.matrix(test.dat)
    
    
    glm.res = glmnet(x=t(train.mat),y=as.factor(train.resp),family="binomial",alpha=1)
    test.res = predict.glmnet(object=glm.res,newx=t(test.mat),y=as.factor(test.resp),type="response")
    lambda = glm.res$lambda
    #lambda = lambdas[,c]
    #lambda = as.vector(lambda)
    #debug
    #print(paste("The length of lambda of current cancer",length(lambda),sep=" "))
    
    prauc.mat = matrix(0,nrow=length(lambda),ncol=2)
    for(i in 1:length(lambda))
    {
      prauc.res = pr.curve(test.res[,i],weights.class0=test.resp.lab,curve=T)
      prauc.mat[i,1] = lambda[i]
      prauc.mat[i,2] = prauc.res$auc.integral
    }
    #colnames(prauc.mat) = c("lambda","prauc")
    tmp.arr.1 = rep(core.cancer[c],length(lambda))
    tmp.arr.2 = rep((cv-3),length(lambda))
    tmp.df = data.frame(tmp.arr.1,tmp.arr.2,prauc.mat)
    colnames(tmp.df) = c("cancer","cv.round","lambda","pr.auc")
    core.res[[cv]] = tmp.df
    #curr.res = rbind(curr.res,tmp.df)
  }
  all.res[[core.cancer[c]]] = core.res
}

#extract lambda for each cancer
lambdas = matrix(0,nrow=100,ncol=4)
for(c in 1:length(core.cancer))
{
  core.res = all.res[[core.cancer[c]]]
  lambda_max = 0
  lambda_max_ix = 0
  for(cv in 4:8)
  {
    if( max(core.res[[cv]]$lambda) > lambda_max )
    {
      lambda_max = max(core.res[[cv]]$lambda)
      lambda_max_ix = cv
    }
  }
  lambdas[,c] = core.res[[lambda_max_ix]]$lambda
}
colnames(lambdas) = core.cancer

#use the unified lambda for each cancer, recalculate
all.res = list()
for(c in 1:length(core.cancer))
{
  core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer[c],]
  core.res = list()
  for(cv in 4:8)
  {
        
    test.pats = as.character(core.info$patient[as.character(core.info[,cv])=="test"])
    test.dat = cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))]
    test.dat = test.dat[,test.pats]
    test.info = core.info[as.character(core.info[,cv])=="test",]
    
    test.resp = as.character(test.info$response[match(colnames(test.dat),as.character(test.info$patient))])
    test.resp.lab = vector(mode="numeric",length=length(test.resp))
    test.resp.lab[which(test.resp=="insensitive")] = 1
    test.resp.lab[which(test.resp=="sensitive")] = 0
    
    train.pats = as.character(core.info$patient[as.character(core.info[,cv])=="train"])
    train.info = core.info[as.character(core.info[,cv])=="train",]
    train.dat = cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))]
    train.dat = train.dat[,train.pats]
    
    #sampling the minority
    #pos.neg.ratio = sum(as.character(core.info$response)=="insensitive")/sum(as.character(core.info$response)=="sensitive")
    pos.neg.ratio = sum(as.character(train.info$response)=="insensitive")/sum(as.character(train.info$response)=="sensitive")
    if( pos.neg.ratio <= 0.6)
    {
      train.minority = as.character(train.info$patient[train.info$response=="insensitive"])
      #sample.num = nrow(train.info) - length(train.minority)
      sample.num = sum(as.character(train.info$response)=="sensitive") - sum(as.character(train.info$response)=="insensitive")
      sample.pats = sample(train.minority,size=sample.num,replace=T)
      col.num = match(sample.pats,colnames(train.dat))
      sample.dat = train.dat[,col.num]
      #train.dat = cbind(train.dat,sample.dat)
      train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
      sample.resp = as.character(train.info$response[match(sample.pats,as.character(train.info$patient))])
      train.resp = c(train.resp,sample.resp)
      train.dat = cbind(train.dat,sample.dat)
      #debug
      train.resp.lab = vector(mode="numeric",length=length(train.resp))
      train.resp.lab[which(train.resp=="insensitive")] = 1
      train.resp.lab[which(train.resp!="insensitive")] = 0
      tmp.char = paste("Current CV round",(cv-3),"of cancer",core.cancer[c],"with",sum(train.resp.lab),"positive samples",sep=" ")
      print(tmp.char)
      
    }
    if( pos.neg.ratio > 0.6)
    {
      train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
      #debug
      train.resp.lab = vector(mode="numeric",length=length(train.resp))
      train.resp.lab[which(train.resp=="insensitive")] = 1
      train.resp.lab[which(train.resp!="insensitive")] = 0
      tmp.char = paste("Current CV round",(cv-3),"of cancer",core.cancer[c],"with",sum(train.resp.lab),"positive samples",sep=" ")
      print(tmp.char)
    }
    
    
    
    train.mat = as.matrix(train.dat)
    test.mat = as.matrix(test.dat)
    
    
    glm.res = glmnet(x=t(train.mat),y=as.factor(train.resp),family="binomial",alpha=1,lambda=lambdas[,c])
    test.res = predict.glmnet(object=glm.res,newx=t(test.mat),y=as.factor(test.resp),type="response")
    #lambda = glm.res$lambda
    lambda = lambdas[,c]
    lambda = as.vector(lambda)
    #debug
    #print(paste("The length of lambda of current cancer",length(lambda),sep=" "))
    
    prauc.mat = matrix(0,nrow=length(lambda),ncol=2)
    for(i in 1:length(lambda))
    {
      prauc.res = pr.curve(test.res[,i],weights.class0=test.resp.lab,curve=T)
      prauc.mat[i,1] = lambda[i]
      prauc.mat[i,2] = prauc.res$auc.integral
    }
    #colnames(prauc.mat) = c("lambda","prauc")
    tmp.arr.1 = rep(core.cancer[c],length(lambda))
    tmp.arr.2 = rep((cv-3),length(lambda))
    tmp.df = data.frame(tmp.arr.1,tmp.arr.2,prauc.mat)
    colnames(tmp.df) = c("cancer","cv.round","lambda","pr.auc")
    core.res[[cv]] = tmp.df
    #curr.res = rbind(curr.res,tmp.df)
  }
  all.res[[core.cancer[c]]] = core.res
}


#start to find the best model in each cancer

best_lambdas = matrix(0,nrow=2,ncol=4)
for(c in 1:length(core.cancer))
{
  curr.res = all.res[[core.cancer[c]]]
  auc_mat = matrix(0,nrow=100,ncol=5)
  auc_mat[,1] = curr.res[[4]]$lambda
  for(cv in 4:8)
  {
    auc_mat[,cv-3] = curr.res[[cv]]$pr.auc
  }
  curr_means = rowMeans(auc_mat)
  curr_ix = which(max(curr_means)==curr_means)
  if(length(curr_ix)>1)
  {
    best_lambdas[1,c] = lambdas[curr_ix[1],c]
    best_lambdas[2,c] = curr_means[curr_ix[1]]
  }
  if(length(curr_ix)==1)
  {
    best_lambdas[1,c] = lambdas[curr_ix,c]
    best_lambdas[2,c] = curr_means[curr_ix]
  }
}
colnames(best_lambdas) = core.cancer




#########################5.5###################################
setwd("C:/Users/zding/workspace/drug_sensitivity/data/omics.drug_centric/cnv")
###input data and information####
cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")

core.cancer = c("LUAD","BLCA","MESO","STAD")
other.cancer = setdiff(as.character(unique(cisplatin.info$cancer)),core.cancer)
library(glmnet)
library(PRROC)

#how to store the data and plot graphs

#train tumor->test tumor -> train model parameters -> test tumor performance
#as list

lambdas = vector(mode="numeric",length=1001)
epsilon = 0.001
lambdas[1001] = 0.5
lambdas[1] = 0.5*epsilon
for(i in 2:1000)
{
  lambdas[i] = 10^0.003 * lambdas[i-1]
}


all.res = list()
for(c in 1:length(core.cancer))
{
  core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer[c],]
  core.res = list()
  for(cv in 4:8)
  {
    
    test.pats = as.character(core.info$patient[as.character(core.info[,cv])=="test"])
    test.dat = cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))]
    test.dat = test.dat[,test.pats]
    test.info = core.info[as.character(core.info[,cv])=="test",]
    
    test.resp = as.character(test.info$response[match(colnames(test.dat),as.character(test.info$patient))])
    test.resp.lab = vector(mode="numeric",length=length(test.resp))
    test.resp.lab[which(test.resp=="insensitive")] = 1
    test.resp.lab[which(test.resp=="sensitive")] = 0
    
    train.pats = as.character(core.info$patient[as.character(core.info[,cv])=="train"])
    train.info = core.info[as.character(core.info[,cv])=="train",]
    train.dat = cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))]
    train.dat = train.dat[,train.pats]
    
    #100 bootstrap samples
    for(i in 1:100)
    {
      
    }
    train.minority = as.character(train.info$patient[train.info$response=="insensitive"])
    #sample.num = nrow(train.info) - length(train.minority)
    sample.num = sum(as.character(train.info$response)=="sensitive") - sum(as.character(train.info$response)=="insensitive")
    sample.pats = sample(train.minority,size=sample.num,replace=T)
    col.num = match(sample.pats,colnames(train.dat))
    sample.dat = train.dat[,col.num]
    #train.dat = cbind(train.dat,sample.dat)
    train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
    sample.resp = as.character(train.info$response[match(sample.pats,as.character(train.info$patient))])
    train.resp = c(train.resp,sample.resp)
    train.dat = cbind(train.dat,sample.dat)
    #debug
    train.resp.lab = vector(mode="numeric",length=length(train.resp))
    train.resp.lab[which(train.resp=="insensitive")] = 1
    train.resp.lab[which(train.resp!="insensitive")] = 0
    tmp.char = paste("Current CV round",(cv-3),"of cancer",core.cancer[c],"with",sum(train.resp.lab),"positive samples",sep=" ")
    print(tmp.char)
   

    
    
    
    train.mat = as.matrix(train.dat)
    test.mat = as.matrix(test.dat)
    
    
    glm.res = glmnet(x=t(train.mat),y=as.factor(train.resp),family="binomial",alpha=1)
    test.res = predict.glmnet(object=glm.res,newx=t(test.mat),y=as.factor(test.resp),type="response")
    lambda = glm.res$lambda
    #lambda = lambdas[,c]
    #lambda = as.vector(lambda)
    #debug
    #print(paste("The length of lambda of current cancer",length(lambda),sep=" "))
    
    prauc.mat = matrix(0,nrow=length(lambda),ncol=2)
    for(i in 1:length(lambda))
    {
      prauc.res = pr.curve(test.res[,i],weights.class0=test.resp.lab,curve=T)
      prauc.mat[i,1] = lambda[i]
      prauc.mat[i,2] = prauc.res$auc.integral
    }
    #colnames(prauc.mat) = c("lambda","prauc")
    tmp.arr.1 = rep(core.cancer[c],length(lambda))
    tmp.arr.2 = rep((cv-3),length(lambda))
    tmp.df = data.frame(tmp.arr.1,tmp.arr.2,prauc.mat)
    colnames(tmp.df) = c("cancer","cv.round","lambda","pr.auc")
    core.res[[cv]] = tmp.df
    #curr.res = rbind(curr.res,tmp.df)
  }
  all.res[[core.cancer[c]]] = core.res
}

#　From github 
