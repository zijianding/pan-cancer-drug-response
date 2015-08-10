####functions###
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
args <- commandArgs(trailingOnly=TRUE)
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")
test_fold = as.numeric(as.character(args[3]))#the current test fold, ranging from 1 to 5
core.cancer = as.character(args[4])

#desktop
# cisplatin.dat = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/Cisplatin.gistic2.gdac_20141206.preprocess.txt",
#                            header=T,row.names=1,sep="\t",quote="")
# cisplatin.info = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt",
#                             sep="\t",header=T,quote="")
# test_fold=1
# core.cancer = "CESC"

#both cluster and desktop
test_fold = test_fold + 3

#for debug
tmp_str = paste("The",test_fold,"test of Cancer",core.cancer,sep=" ")
print(tmp_str)
#

#libraries
library(randomForest)
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
models = c(seq(0.01,0.5,by=0.01),1) #features, "1" means all features(bagging)
folds = 5 # 5-fold CV for best model selection
BS = 100 # 100 bootstrap samples one model
score.mat = matrix(0,nrow=length(models),ncol=length(train.pats)) 
colnames(score.mat) = train.pats
cv.mat = partition_data(train.pats,5)

for(fold in 1:folds)
{
  #print messegae
  tmp_str = paste("Current iteration:",fold,"of 5-fold cross validation",sep=" ")
  print(tmp_str)
  ###prepare data###
  #find current validation sets
  curr.validation_pats = rownames(cv.mat)[cv.mat[,fold]=="validation"]
  curr.validation_resp = as.character(train.info$response[match(curr.validation_pats,as.character(train.info$patient))])
  curr.validation_dat =  t( cisplatin.dat[,match(curr.validation_pats,as.character(colnames(cisplatin.dat)))] ) 
  #find current train sets
  curr.train_pats = rownames(cv.mat)[cv.mat[,fold]=="train"]
  curr.train_resp = as.character(train.info$response[match(curr.train_pats,as.character(train.info$patient))])
  curr.train_dat = t( cisplatin.dat[,match(curr.train_pats,colnames(cisplatin.dat))] ) 
  #curr.train_dat = data.frame(curr.train_dat,curr.train_resp)
  
  all.score = array( rep(0,length(curr.validation_resp)*length(models)*BS),
                     c(length(curr.validation_resp),length(models),BS) )
  for(bs in 1:BS)
  {
    #sample the minority class
    sample.num = floor(length(curr.train_resp)/2)
    random.sample = sample(curr.train_pats[curr.train_resp=="insensitive"],replace=T,
                           size=length(curr.train_resp)-sum(curr.train_resp=="insensitive"))
    random.resp = curr.train_resp[match(random.sample,curr.train_pats)]
    random.dat = t( cisplatin.dat[,match(random.sample,as.character(colnames(cisplatin.dat)))])
    train_dat = rbind(curr.train_dat,random.dat)
    train_resp = c(curr.train_resp,random.resp)
    
    #start random forest model
    for(i in 1:length(models))
    {
      rf.fit = randomForest(x=train_dat,y=as.factor(train_resp),ntree=BS,replace=T,
                            strata=as.factor(train_resp),mtry=floor(models[i]*ncol(train_dat)),
                            sampsize=c(sample.num,sample.num) )
      pred = predict(object=rf.fit,newdata=curr.validation_dat,type="prob")
      
      ix = match("insensitive",colnames(pred))
      all.score[,i,bs] = pred[,ix]
    }
  }
  
  #average the output
  test.ix = match(curr.validation_pats,colnames(score.mat))
  for(i in 1:length(test.ix))
  {
    curr_score = all.score[i,,]
    curr_score_ave = rowMeans(curr_score)
    score.mat[,test.ix[i]] = as.vector(curr_score_ave)
  }
  
}

####best model performance of CV####
#calculate PR curve and ROC curve
auc.best = 0
auc.best_ix = 0
#change to parallel
cl = makeCluster(no_cores)
registerDoParallel(cl)
auc.score<-foreach(i=1:length(models),.combine='c') %dopar%
{
  library(PRROC)
  prauc = pr.curve(score.mat[i,],weights.class0=train.resp.lab,curve=T)
  return(prauc$auc.integral)
}
stopImplicitCluster()
stopCluster(cl)
auc.best = max(auc.score)
auc.best_ix = which.max(auc.score)
pr_best = as.vector(models[auc.best_ix])
prauc.best = pr.curve(score.mat[auc.best_ix,],weights.class0=train.resp.lab,curve=T)
tmp_str = paste("Single cancer",core.cancer,"by Random Forest(CV)",sep=" ")
plot(prauc.best,main=tmp_str);
legend( "topright",legend=paste("Number of features:",pr_best,"%",sep="") )

#classification score
auc.best = 0
auc.best_ix = 1
cl = makeCluster(no_cores)
registerDoParallel(cl)
auc.score <-foreach(i=1:length(models),.combine='c') %dopar%
{
  library(PRROC)
  roc = roc.curve(score.mat[i,],weights.class0=train.resp.lab,curve=T)
  return(roc$auc)
}
stopImplicitCluster()
stopCluster(cl)
auc.best = max(auc.score)
auc.best_ix = which.max(auc.score)
roc_best = models[auc.best_ix]
roc.best = roc.curve(score.mat[auc.best_ix,],weights.class0=train.resp.lab,curve=T)
tmp_str = paste("Single cancer",core.cancer,"by random forest(CV)", sep=" ")
plot(roc.best,main=tmp_str);
legend( "bottomright",legend=paste("Number of features",roc_best,"%",sep="") )


####test performance of best models####
#best PR and ROC curve model together
feature_gini.pr = matrix(0,nrow=nrow(test.dat),ncol=BS)#temporarily no idea
feature_gini.roc = matrix(0,nrow=nrow(test.dat),ncol=BS)#temporarily no idea
#1 best PR model, 100 BS samples,therefore 100 scores for each observation
array_score.pr = matrix(0,nrow=length(test.pats),ncol=100)
array_score.roc = matrix(0,nrow=length(test.pats),ncol=100) 
for(bs in 1:BS)
{
  #make even the class
  sample.pats = sample( train.pats[train.resp=="insensitive"],replace=T,
                       size=length(train.resp)-sum(train.resp=="insensitive") )
  sample.resp = rep("insensitive",length(sample.pats))
  sample.num = floor(length(sample.pats)/2)
  sample.dat = t( cisplatin.dat[,match(sample.pats,as.character(colnames(cisplatin.dat)))] )
  
  curr.train_pats = c(train.pats,sample.pats)
  curr.train_resp = c(train.resp,sample.resp)
  curr.train_dat = rbind(t(train.dat),sample.dat)
  #fit pr best model to all training data and test
  rf_fit.pr = randomForest(x=curr.train_dat,y=as.factor(curr.train_resp),ntree=BS,
                        strata=as.factor(train_resp),mtry=floor(pr_best*ncol(curr.train_dat)),
                        sampsize=c(sample.num,sample.num),importance=T )
  
  feature_gini.pr[,bs] = rf_fit.pr$importance[,4] #not totally understand

  pred.pr = predict(object=rf_fit.pr,newdata=t(test.dat),type="prob")
  ix = match( "insensitive",colnames(pred.pr) )
  array_score.pr[,bs] = pred.pr[,ix]
  
  #fit roc best model all training data and test
  #fit pr best model to all training data and test
  rf_fit.roc = randomForest(x=curr.train_dat,y=as.factor(curr.train_resp),ntree=BS,
                        strata=as.factor(train_resp),mtry=floor(roc_best*ncol(curr.train_dat)),
                        sampsize=c(sample.num,sample.num),importance=T )
  feature_gini.roc[,bs] = rf_fit.roc$importance[,4] #not totally understand
  pred.roc = predict(object=rf_fit.roc,newdata=t(test.dat),type="prob")
  ix = match("insensitive",colnames(pred.roc))
  array_score.roc[,bs] = pred.roc[,ix]
  
}
classification_score.pr = rowMeans(array_score.pr)
classification_score.roc = rowMeans(array_score.roc)

gini.pr = rowMeans(feature_gini.pr)
gini.roc = rowMeans(feature_gini.roc)

#plot pr and roc curve, and corresponding AUC
#PR curve
tmp_str = paste("Single cancer",core.cancer,"by Random Forest(TEST)",sep=" ")
prauc.best = pr.curve(classification_score.pr,weights.class0=test.resp.lab,curve=T)
plot(prauc.best,main=tmp_str);
legend( "topright",legend=paste("Number of features:",pr_best,"%",sep="") )

#ROC curve
roc.best = roc.curve(classification_score.roc,weights.class0=test.resp.lab,curve=T)
tmp_str = paste("Single cancer",core.cancer,"by random forest(TEST)", sep=" ")
plot(roc.best,main=tmp_str);
legend( "bottomright",legend=paste("Number of features",roc_best[2],"%",sep="") )


####output characteristics of best model####
#output the classification score
classify_score = cbind(classification_score.pr,classification_score.roc)
colnames(classify_score) = c("PR","ROC")
rownames(classify_score) = test.pats
tmp_str = paste("classify_score.randomForest",core.cancer,"test_",test_fold-3,"txt",sep=".")
write.table(classify_score,file=tmp_str,row.names=T,col.names=T,sep="\t",quote=F)

#output the feature importance
gini = cbind(gini.pr,gini.roc)
colnames(gini) = c("PR","ROC")
rowneams(gini) = rownames(train.dat)
tmp_str = paste("gini.randomForest",core.cancer,"test_",test_fold-3,"txt",sep=".")
write.table(gini,file=tmp_str,row.names=T,col.names=T,sep="\t",quote=F)






