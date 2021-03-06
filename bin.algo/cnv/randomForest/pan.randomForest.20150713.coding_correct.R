####for CNV data###
###system parameters###
#options(digits = 15)

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

#for debug
# train_dat=train.dat; test_dat = test.dat;info=cisplatin.info 
# type="regression"; parallel=T
#
test_gene <- function(train.dat, test.dat, sig_gene=50, p_thresh=0.05,type="ttest")
{
  #data are generated by
  if(type=="ttest")
  {
    pvalue = vector(length=ncol(train.dat)-1,mode="numeric")
    ix.1 = which(train.dat$Class=="Class1")
    ix.2 = which(train.dat$Class=="Class2")
    for( i in 1:(ncol(train.dat)-1) )
    {
      dat.1 = train.dat[ix.1,i]
      dat.2 = train.dat[ix.2,i]
      test_fit = t.test(dat.1,dat.2,na.action="na.omit")
      pvalue[i] = test_fit$p.value
    }
    ix = which(pvalue<=0.05)
    ix = c(ix,ncol(train.dat))
    return( list(train.dat[,ix],test.dat[,ix],length(ix)-1) )
  }
}

gene_selection <- function(data,freq=0.8)
{
  #data: row as BS and column as genes
  #new_data = data
  #new_data[new_data!=0] = 1
  
  #direction
  data[data>0] = 1
  data[data<0] = -1
  #times of larger/smaller than 0
  up = colSums(data>0)
  down = colSums(data<0)
  #the largest times and freq
  side_max = abs(up-down)
  side_freq = side_max/nrow(data)
  names(side_freq) = colnames(data)
  
  
  #feature_freq = colSums(new_data)
  #feature_freq = feature_freq/nrow(data)
  
  
  ix = which( side_freq>=freq )
  
  genes = colnames(data)[ix]
  
  return(list(genes,ix,side_freq))
  
}

ensemble_roc <- function(score_mat,class,target_class,step=0.0001)
{
  #score_mat: row as one model, column as patients
  #class: responses
  #return value: tpr and fpr under each cut
  #cut: seq(0,1,by=0.01)
  #currently target class is treated as negative classs
  #which is wrong, later we will corrected this
  # for debug
  #     score_mat = t(test_score)
  #     class = test.resp
  #     target_class = "sensitive"
  #     step=0.0001
  
  #
  
  
  #here is the problem
  class_1 = target_class # large score
  class_2 = setdiff( unique(class), target_class )
  #   if( sum(class==class_1) > sum(class==class_2) )
  #   {
  #     class_p = class_2
  #     class_n = class_1
  #   }
  #   if( sum(class==class_1) <= sum(class==class_2) )
  #   {
  #     class_p = class_1
  #     class_n = class_2
  #   }
  class_p = class_2 #small score
  class_n = class_1 #large score
  
  
  
  cut_vec = seq(0,1+step,by=step)
  cut_vec[length(cut_vec)] = Inf
  tpr_vec = vector(length=length(cut_vec),mode="numeric")
  fpr_vec = vector(length=length(cut_vec),mode="numeric")
  
  for( i in 1:length(cut_vec) )
  {
    curr_mat = matrix( NA,nrow=nrow(score_mat),ncol=ncol(score_mat) )
    curr_mat[ score_mat>=cut_vec[i] ] = 1 #large score
    curr_mat[ score_mat< cut_vec[i] ] = 0 #small score
    
    curr_score = colMeans( curr_mat )
    
    curr_class = vector(length=length(curr_score),mode="character")
    curr_class[curr_score>0.5] = class_1
    curr_class[curr_score<=0.5] = class_2
    
    #mytable = table( curr_class,class )
    mytable = matrix(0,nrow=2,ncol=2)
    dimnames(mytable) = list( c(class_p,class_n), c(class_p,class_n) )
    
    mytable[1,1] = sum(curr_class[class==class_p]==class_p)
    mytable[2,1] = sum(curr_class[class==class_p]==class_n)
    mytable[1,2] = sum(curr_class[class==class_n]==class_p)
    mytable[2,2] = sum(curr_class[class==class_n]==class_n)
    
    
    if( sum(mytable[,1])!=0 )
    {
      tpr_vec[i] = mytable[1,1]/sum(mytable[,1])
    }
    if( sum(mytable[,1])==0 )
    {
      tpr_vec[i] = NA 
    }
    if( sum(mytable[,2])!=0 )
    {
      fpr_vec[i] = mytable[1,2]/sum(mytable[,2])
    }
    if( sum(mytable[,2])==0 )
    {
      fpr_vec[i] = NA
    }
    
  }
  
  roc_mat = cbind(tpr_vec,fpr_vec)
  colnames(roc_mat) = c("tpr","fpr")
  return(roc_mat)
  
}


####load data####
#libraries
library(glmnet)
library(doParallel)
library(foreach)
library(pracma)
library(caret)
no_cores = detectCores()


##find data##
#artifical test data#
library(caret)
set.seed(1)
#dat <- twoClassSim(1000, intercept = -50,linearVars=1000,noiseVars=50,corrVars=10)
dat <- twoClassSim(1000, intercept = -10)
table(dat$Class)
datpart = createDataPartition(y=dat$Class, times=1, p=0.5, list=F)
train.dat = dat[datpart[,1],]
test.dat = dat[-datpart[,1],]
#test.dat <- twoClassSim(5000, intercept = -50,linearVars=1000,noiseVars=50,corrVars=10)


##select differential genes across cancer types##
list_tmp = test_gene(train.dat, test.dat)
feature_num = list_tmp[[3]]
print(feature_num)
train.dat = list_tmp[[1]]
test.dat = list_tmp[[2]]



###estimated robust features###
#basic parameters#
BS = 100

#bootstrap samples#
#train patients response
test.resp = as.character(test.dat$Class)
train.resp = as.character(train.dat$Class)
min.class = "Class2"
sampsize = ceiling(length(train.resp)/2)
min.ix = which(train.resp==min.class)
samptime = floor(sampsize/length(min.ix))
for(i in 1:samptime)
{
  train.dat = rbind( train.dat, train.dat[min.ix,])
  train.resp = c( train.resp,rep(min.class,times=length(min.ix)) )
}

##estimates recurrent features##
#model training#

ctrl <- trainControl(method = "cv",
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
set.seed(2)
cl = makeCluster(no_cores)
registerDoParallel(cl)
rf <- train(x=train.dat[,-ncol(train.dat)], 
            y=as.factor(train.resp),
            method = "rf",
            ntree = BS,
            tuneLength = 10,
            metric = "ROC",
            trControl = ctrl,
            ## Tell randomForest to sample by strata. Here, 
            ## that means within each class
            strata = as.factor(train.resp),
            ## Now specify that the number of samples selected
            ## within each class should be the same
            sampsize = rep(sampsize, 2),
            importance=T)
                       
stopImplicitCluster()
stopCluster(cl)

#find best features
fm = rf$finalModel
feature_num = rf$bestTune[[1]]
ix = sort(fm$importance[,ncol(fm$importance)],decreasing=T,index.return=T)$ix[1:feature_num]
tmp_str = paste( "Select recurrent features",feature_num )
print(tmp_str)

#refine data by recurrent features
train.dat = train.dat[,ix]
test.dat = test.dat[,ix]



##refit models and test##
bag_tree = randomForest(x=train.dat,y=as.factor(train.resp),
                        ntree=BS, mtry=feature_num,
                        strata=as.factor(train.resp),
                        sampsize=rep(sampsize,2))

pred = predict(object=bag_tree,newdata=test.dat,type="prob")


###output results###
#roc
ROC <- roc(response = test.resp, 
           predictor = pred[,2])
plot(ROC,xlab="Specificity",ylab="Sensitivity")
title("Pan-cancer performance",paste("AUC =",round(ROC$auc,digits=2)))
dev.off()

####THE END####

            
