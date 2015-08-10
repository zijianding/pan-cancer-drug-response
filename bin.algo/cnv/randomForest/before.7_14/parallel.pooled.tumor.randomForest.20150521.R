####functions###
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

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}
####load data###
#cluster
args <- commandArgs(trailingOnly=TRUE)
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")
test_fold = as.numeric(as.character(args[3]))#the current test fold, ranging from 1 to 10
output_folder = args[4] #no "/" at the end

#desktop
cisplatin.dat = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/Cisplatin.gistic2.gdac_20141206.preprocess.txt",
                           header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt",
                            sep="\t",header=T,quote="")
test_fold=1
output_folder = "."

#both cluster and desktop
test_fold = test_fold + 3

#libraries
library(randomForest)
library(PRROC)
library(doParallel)
library(foreach)
no_cores = detectCores()

###preprocess data###
core.info = cisplatin.info
core.cancer = c("CESC","LUAD", "BLCA")
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
cv.mat = partition_data(train.pats,folds)

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
  
  
  
  list.bs_mat = bootstrap_sample(curr.train_pats,curr.train_resp,BS,method="supply")
  bs_mat.pats = list.bs_mat[[1]]
  bs_mat.resp = list.bs_mat[[2]]
  
  #a parallel version
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  curr_res <- foreach( i=rep(1:length(models),each=BS),
                       bs=rep(1:BS,times=length(models)) ) %dopar%
  {

	  curr.train_dat = t( cisplatin.dat[,match(bs_mat.pats[,bs],colnames(cisplatin.dat))] )
    sample.num = floor(length(curr.train_resp)/2)
	
	  #start random forest model
	  library(randomForest)
    rf.fit = randomForest(x=curr.train_dat,y=as.factor(bs_mat.resp[,bs]),ntree=BS,replace=T,
	                        strata=as.factor(bs_mat.resp[,bs]), mtry=floor(models[i]*ncol(curr.train_dat)),
						              sampsize=c(sample.num,sample.num))
	  pred = predict(object=rf.fit,newdata=curr.validation_dat,type="prob")
	  ix = match("insensitive",colnames(pred))
	  return(pred[,ix])
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  #   all.score = array( rep(0,length(curr.validation_resp)*length(models)*BS),
  #                      c(length(curr.validation_resp),length(models),BS) )
  all.score = array(data=NA,dim=c(length(curr.validation_resp),length(models),BS) )
  
  #position the output
  for(i in 1:length(curr_res) )
  {
    model_ix = ceiling(i/BS)
    bs = i%%BS
    if(bs==0)
    {
      bs=BS
      all.score[,model_ix,bs] = curr_res[[i]]
    }
    if(bs!=0)
    {
      all.score[,model_ix,bs] = curr_res[[i]]
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

# #for debug
tmp_str = paste(output_folder,"/","Pan(pooled).test_",test_fold-3,".randomForest.score_mat(CV).txt",sep="")
write.table(score.mat,tmp_str,col.names=T,row.names=T,quote=F,sep="\t")
# score.mat = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/randomForest/pool/Pan(pooled).test_1.randomForest.score_mat(CV).txt",
#                        header=T,sep="\t",quote="",row.names=1)
# score.mat = as.matrix( score.mat[,match(train.pats,colnames(score.mat))] )
# #

####best model performance of CV####
#calculate PR/ROC curve
auc_score.mat = matrix(0,nrow=length(models),ncol=(2+2*length(core.cancer)))
colnames(auc_score.mat) = c("pr",core.cancer,"roc",core.cancer)
for(i in 1:nrow(score.mat))
{
  ##pan-cancer##
  #pr auc 
  prauc = pr.curve(score.mat[i,],weights.class0=train.resp.lab,curve=T)
  auc_score.mat[i,1] = prauc$auc.integral
  #roc auc
  roc = roc.curve(score.mat[i,],weights.class0=train.resp.lab,curve=T)
  auc_score.mat[i,2+length(core.cancer)] = roc$auc
  
  ##on core cancer##
  for(j in 1:length(core.cancer))
  {
    #find patients
    curr.pats = train.pats[core.info$cancer[match(train.pats,as.character(core.info$patient))] == core.cancer[j]]
    curr.ix = match( curr.pats,colnames(score.mat) )
    curr.score = as.vector(score.mat[i,curr.ix])
    curr.resp_lab = train.resp.lab[curr.ix]
    #pr
    prauc = pr.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
    auc_score.mat[i,j+1] = prauc$auc.integral
    #roc
    roc = roc.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
    auc_score.mat[i,2+length(core.cancer)+j] = roc$auc
  }
}
#best value and paramters
auc_best.mat = matrix(0,nrow=3,ncol=2+2*length(core.cancer) )
dimnames(auc_best.mat) = list(c("auc","param","param.ix"),
                              c("pr",core.cancer,"roc",core.cancer))
for( j in 1:ncol(auc_score.mat))
{
  curr.best = max(auc_score.mat[,j])
  curr.best_ix = which.max(auc_score.mat[,j])
  curr.param = models[curr.best_ix]
  
  auc_best.mat[1,j] = curr.best 
  auc_best.mat[2,j] = curr.param
  auc_best.mat[3,j] = curr.best_ix
}
#store CV performance
cv_score.pr_roc = matrix(0,nrow=(2+2*length(core.cancer)),ncol=length(train.pats))
dimnames(cv_score.pr_roc) = list(c("pr",core.cancer,"roc",core.cancer),train.pats)
best_cv.mat = matrix(0, nrow=2+2*length(core.cancer),ncol=ncol(score.mat))
dimnames(best_cv.mat) = list(c("pr",core.cancer,"roc",core.cancer),colnames(score.mat))
for( j in 1:nrow(best_cv.mat))
{
  cv_score.pr_roc[j,] = score.mat[auc_best.mat[3,j],]
}
#draw the CV performance
file_name = paste(output_folder,"/","performance.randomForest.Pan(pooled).test_",test_fold-3,".20150526.pdf",sep="")
pdf(file=file_name)
##pan-cancer
#pr curve, as a whole
prauc = pr.curve(cv_score.pr_roc[1,],weights.class0=train.resp.lab,curve=T)
tmp_str = "Pan-cancer(PR) by Random Forest(CV)"
plot(prauc,main=tmp_str);
legend( "topright",legend=paste("Number of features:",auc_best.mat[2,1],sep="") )
#each cancer
for(i in 1:length(core.cancer))
{
  #find patients
  curr.pats = train.pats[core.info$cancer[match(train.pats,as.character(core.info$patient))] == core.cancer[i]]
  curr.ix = match( curr.pats,colnames(score.mat) )
  curr.score = as.vector(cv_score.pr_roc[1,curr.ix])
  curr.resp_lab = train.resp.lab[curr.ix]
  #best on single
  tmp_str = paste("Pan-cancer(PR) performance(CV) by Random Forest on", core.cancer[i],sep=" ")
  prauc = pr.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
  plot(prauc,main=tmp_str);
  legend( "topright",legend=paste("Number of features:",auc_best.mat[2,1],sep="") )
  #single best
  curr.score = as.vector(cv_score.pr_roc[1+i,curr.ix])
  tmp_str = paste("Pan-cancer(PR) performance(CV) by Random Forest chosen by", core.cancer[i],sep=" ")
  prauc = pr.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
  plot(prauc,main=tmp_str);
  legend( "topright",legend=paste("Number of features:",auc_best.mat[2,1+i],sep="") )
}
#roc curve, as a whole
roc.best = roc.curve(cv_score.pr_roc[2+length(core.cancer),],weights.class0=train.resp.lab,curve=T)
tmp_str = "Pan-cancer(ROC) by Random Forest(CV)"
plot(roc.best,main=tmp_str);
legend( "bottomright",legend=paste("Number of features:",auc_best.mat[2,2+length(core.cancer)],sep="") )
#each cancer
for(i in 1:length(core.cancer))
{
  #find patients
  curr.pats = train.pats[core.info$cancer[match(train.pats,as.character(core.info$patient))] == core.cancer[i]]
  curr.ix = match( curr.pats,colnames(score.mat) )
  curr.score = as.vector(cv_score.pr_roc[2+length(core.cancer),curr.ix])
  curr.resp_lab = train.resp.lab[curr.ix]
  #best on single
  tmp_str = paste("Pan-cancer(PR) performance(CV) by Random Forest on", core.cancer[i],sep=" ")
  prauc = roc.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
  plot(prauc,main=tmp_str);
  legend( "topright",legend=paste("Number of features:",auc_best.mat[2,2+length(core.cancer)],sep="") )
  #single best
  curr.score = as.vector( cv_score.pr_roc[2+length(core.cancer)+i,curr.ix] )
  tmp_str = paste("Pan-cancer(PR) performance(CV) by Random Forest chosen by", core.cancer[i],sep=" ")
  prauc = roc.curve( curr.score,weights.class0=curr.resp_lab,curve=T )
  plot(prauc,main=tmp_str);
  legend( "topright",legend=paste("Number of features:",auc_best.mat[2,2+length(core.cancer)+i],sep="") )
}



####test performance of best models####
#best PR&ROC curve models together, #(2+2*length(core.cancer)) models
feature_gini.arr = array( data=NA,dim=c(nrow(test.dat),2*length(core.cancer)+2,BS) )#temporarily no idea right or wrong

#1 best PR&ROC model, 100 BS samples,therefore 100 scores for each observation
score.arr = array( data=NA,c(length(test.pats),2*length(core.cancer)+2,BS) )

#bootstrap the minority class
list.bs_mat = bootstrap_sample(train.pats,train.resp,BS,method="supply")
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]

for(i in 1:ncol(auc_best.mat))
{
  curr.mtry = auc_best.mat[2,i]
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  curr_res <- foreach(bs=1:BS,.combine='comb', .multicombine=TRUE,
                      .init=list(list(), list()) ) %dopar%
  {
    sample.num = floor(length(train.pats)/2)
    curr.train_pats = bs_mat.pats[,bs]
    curr.train_resp = bs_mat.resp[,bs]
    curr.train_dat = t(cisplatin.dat[,match(curr.train_pats,as.character(colnames(cisplatin.dat)))])
    #fit pr best model to all training data and test
    library(randomForest)
    rf_fit.pr = randomForest(x=curr.train_dat,y=as.factor(curr.train_resp),ntree=BS,
                             strata=as.factor(curr.train_resp),mtry=floor(curr.mtry*ncol(curr.train_dat)),
                             sampsize=c(sample.num,sample.num),importance=T )
    
    #feature_gini.arr[,i,bs] = rf_fit.pr$importance[,4] #not totally understand
    curr_fg = rf_fit.pr$importance[,4]
    pred.pr = predict(object=rf_fit.pr,newdata=t(test.dat),type="prob")
    ix = match( "insensitive",colnames(pred.pr) )
    #score.arr[,i,bs] = pred.pr[,ix]
    curr_score = pred.pr[,ix]
    
    return(list(curr_fg,curr_score))
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  #value assignment
  feature_gini = curr_res[[1]]
  score = curr_res[[2]]
  for(bs in 1:BS)
  {
    feature_gini.arr[,i,bs] = feature_gini[[bs]]
    score.arr[,i,bs] = score[[bs]]
  }
  
}

#average classification score
ave_score.arr = matrix(0,nrow=length(test.pats),ncol=ncol(auc_best.mat))
ave_feature_gini.arr = matrix(0,nrow=nrow(test.dat),ncol=ncol(auc_best.mat))
for(i in 1:ncol(auc_best.mat))
{
  ave_score.arr[,i] = rowMeans(score.arr[,i,])
  ave_feature_gini.arr = rowMeans(feature_gini.arr[,i,])
}
dimnames(ave_score.arr) = list( test.pats,colnames(auc_best.mat) )
names(ave_feature_gini.arr) = rownames(test.dat)

#plot pr and roc curve, and corresponding AUC
for(i in 1:ncol(auc_best.mat))
{
  if(colnames(auc_best.mat)[i] == "pr" )
  {
    #pan-cancer
    tmp_str = "Pan-cancer(PR) by Random Forest(TEST)"
    prauc.best = pr.curve(ave_score.arr[,i],weights.class0=test.resp.lab,curve=T)
    plot(prauc.best,main=tmp_str);
    legend( "topright",legend=paste("Number of features:",auc_best.mat[2,i],sep="") )
    
  }
  if( colnames(auc_best.mat)[i] == "roc" )
  {
    #pan-cancer
    roc.best = roc.curve(ave_score.arr[,i],weights.class0=test.resp.lab,curve=T)
    tmp_str = "Pan-cancer(ROC) by Random Forest(TEST)"
    plot(roc.best,main=tmp_str);
    legend( "bottomright",legend=paste("Number of features:",auc_best.mat[2,i],sep="") )
    
  }
  if(colnames(auc_best.mat)[i] %in% core.cancer )
  {
    #corresponding patients
    curr.pats = test.pats[core.info$cancer[match(test.pats,as.character(core.info$patient))] == colnames(auc_best.mat)[i] ]
    curr.ix = match(curr.pats,rownames(ave_score.arr))
    curr.score = as.vector(ave_score.arr[curr.ix,i])
    curr.resp_lab = test.resp.lab[curr.ix]
    # PR or ROC
    if(i<=(1+length(core.cancer)))
    {
      #PR
      tmp_str = paste("Pan-cancer(PR) by Random Forest(TEST) on",colnames(auc_best.mat)[i],sep=" ")
      prauc.best = pr.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
      plot(prauc.best,main=tmp_str);
      legend( "topright",legend=paste("Number of features:",auc_best.mat[2,i],sep="") )
    }
    if(i>(1+length(core.cancer)))
    {
      #ROC
      roc.best = roc.curve(curr.score,weights.class0=curr.resp_lab,curve=T)
      tmp_str = paste("Pan-cancer(ROC) by Random Forest(TEST) on",colnames(auc_best.mat)[i],sep=" ")
      plot(roc.best,main=tmp_str);
      legend( "bottomright",legend=paste("Number of features:",auc_best.mat[2,i],sep="") )
    }
    
  }
}






####output characteristics of best model####
#output the best models
tmp_str = paste(output_folder,"/pan(pool).param.randomForest.test_",test_fold-3,".txt",sep="")
write.table(auc_best.mat,file=tmp_str,row.names=T,col.names=T,sep="\t",quote=F)

#output the CV performance
tmp_str = paste(output_folder,"/","pan(pool).cv_score.randomForest.test_",test_fold-3,".txt",sep="")
write.table(cv_score.pr_roc,file=tmp_str,row.names=T,col.names=T,sep="\t",quote=F)

#output test performance
tmp_str = paste(output_folder,"/","pan(pool).test_score.randomForest.test_",test_fold-3,".txt",sep="")
write.table(ave_score.arr,file=tmp_str,row.names=T,col.names=T,sep="\t",quote=F)

#output the feature importance
tmp_str = paste(output_folder,"/","pan(pooled).gini.randomForest.test_",test_fold-3,".txt",sep="")
write.table(ave_feature_gini.arr,file=tmp_str,row.names=T,col.names=T,sep="\t",quote=F)


dev.off()



