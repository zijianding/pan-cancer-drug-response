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
test_gene <- function(train_dat, test_dat, info, sig_gene = 50,
                      p_thresh=0.05,q_thresh=0.05, type="ttest",parallel=F)
{
  #input data: row as genes and col as samples
  #output data:row as genes and col as samples
  #output significant threshold: q_thresh and p_thresh
  #type: default is t test; 
  #wilcox rank sum test, and logistic regression are also provided: "wilcox", "regression"
  #sig_gene, when no significant genes according to p/q value threshold
  
  
  if( type=="ttest" )
  {
    #perform on train data
    train_pats = colnames(train_dat)
    curr_info = info[match(train_pats,as.character(info$patient)),]
    sen_ix = which(as.character(curr_info$response)=="sensitive")
    insen_ix = which(as.character(curr_info$response)=="insensitive")
    
    #perform test
    p_values = vector( length=nrow(train_dat),mode="numeric" )
    if( parallel==F)
    {
      for( i in 1:nrow(train_dat) )
      {
        sen_dat = as.numeric(train_dat[i,sen_ix])
        insen_dat = as.numeric(train_dat[i,insen_ix])
        test_fit = t.test(sen_dat,insen_dat,na.action="na.omit")
        p_values[i] = test_fit$p.value
      }
    }
    if( parallel==T)
    {
      p_values <- foreach( i=1:nrow(train_dat),
                           .combine='c') %dopar%{
                             
                             sen_dat = as.numeric(train_dat[i,sen_ix])
                             insen_dat = as.numeric(train_dat[i,insen_ix])
                             test_fit = t.test(sen_dat,insen_dat,na.action="na.omit")
                             return(test_fit$p.value)
                           }
    }
    
    
    #multiple correction
    q_values = p.adjust(p_values,method="fdr",n=length(p_values))
    
    #determine q value threshold or p value threshold
    q_use = TRUE
    while( sum(q_values<q_thresh) <= sig_gene )
    {
      q_thresh = q_thresh + 0.05
      if( q_thresh > 0.25)
      {
        q_use = FALSE
        break
      }
    }
    q_use = FALSE #use p value at last
    p_use = TRUE
    if( q_use == FALSE )
    {
      while( sum(p_values<p_thresh) <= sig_gene  )
      {
        p_thresh = p_thresh + 0.05
        if(p_thresh > 0.15)
        {
          p_use = FALSE
          break
        }
      }
    }
    
    #identify significant genes
    if( q_use == TRUE )
    {
      sig_ix = q_values<=q_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"q_thresh",q_thresh))
    }
    if( (q_use==F) && (p_use==T) )
    {
      sig_ix = p_values <= p_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"p_thresh",p_thresh))
    }
    if( (q_use==F) && (p_use==F) )
    {
      sig_ix = sort(p_values,decreasing=F,index.return=T)$ix
      sig_ix = sig_ix[1:sig_gene]
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"sig_gene",sig_gene))
    }
    
  }
  if( type=="wilcox" )
  {
    #train data
    train_pats = colnames(train_dat)
    curr_info = info[match(train_pats,as.character(info$patient)),]
    sen_ix = which(as.character(curr_info$response)=="sensitive")
    insen_ix = which(as.character(curr_info$response)=="insensitive")
    
    #perform test
    p_values = vector( length=nrow(train_dat),mode="numeric" )
    if( parallel==F )
    {
      for( i in 1:nrow(train_dat) )
      {
        sen_dat = as.numeric(as.character(train_dat[i,sen_ix]))
        insen_dat = as.numeric(as.character(train_dat[i,insen_ix]))
        test_fit = wilcox.test(sen_dat,insen_dat,na.action="na.omit")
        p_values[i] = test_fit$p.value
      }
    }
    if( parallel==T)
    {
      p_values <- foreach( i=1:nrow(train.dat), .combine='c' ) %dopar%{
        sen_dat = as.numeric(as.character(train_dat[i,sen_ix]))
        insen_dat = as.numeric(as.character(train_dat[i,insen_ix]))
        test_fit = wilcox.test(sen_dat,insen_dat,na.action="na.omit")
        return(test_fit$p.value)
      }
    }
    
    #multiple correction
    q_values = p.adjust(p_values,method="fdr",n=length(p_values))
    
    #determine q value threshold or p value threshold
    q_use = TRUE
    while( sum(q_values<q_thresh) <= sig_gene )
    {
      q_thresh = q_thresh + 0.05
      if( q_thresh > 0.25)
      {
        q_use = FALSE
        break
      }
    }
    q_use = FALSE #use p value at last
    p_use = TRUE
    if( q_use == FALSE )
    {
      while( sum(p_values<p_thresh) <= sig_gene  )
      {
        p_thresh = p_thresh + 0.05
        if(p_thresh > 0.15)
        {
          p_use = FALSE
          break
        }
      }
    }
    
    #identify significant genes
    if( q_use == TRUE )
    {
      sig_ix = q_values<=q_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"q_thresh",q_thresh))
    }
    if( (q_use==F) && (p_use==T) )
    {
      sig_ix = p_values <= p_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"p_thresh",p_thresh))
    }
    if( (q_use==F) && (p_use==F) )
    {
      sig_ix = sort(p_values,decreasing=F,index.return=T)$ix
      sig_ix = sig_ix[1:sig_gene]
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"sig_gene",sig_gene))
    }
    
  }
  if( type =="regression")
  {
    #train data
    train_pats = colnames(train_dat)
    cancers = as.character(info$cancer[match(train_pats,as.character(info$patient))])
    cancers = as.factor(cancers)
    responses = as.character(info$response[match(train_pats,as.character(info$patient))])
    responses = as.factor(responses)
    
    #perform test
    p_values = vector( length=nrow(train_dat),mode="numeric" )
    if( parallel==F )
    {
      for(i in 1:nrow(train_dat) )
      {
        fit = glm(responses~as.numeric(as.character(train_dat[i,]))+cancers,family=binomial)
        sum_fit = summary(fit)
        coeff = sum_fit$coefficients
        p_values[i] = coeff[2,ncol(coeff)]
      }
    }
    if( parallel==T )
    {
      p_values <- foreach(i=1:nrow(train_dat), .combine='c') %dopar%{
        fit = glm(responses~as.numeric(as.character(train_dat[i,]))+cancers,family=binomial)
        sum_fit = summary(fit)
        coeff = sum_fit$coefficients
        return( coeff[2,ncol(coeff)] )
      }
      
    }
    
    #multiple correction
    q_values = p.adjust(p_values,method="fdr",n=length(p_values))
    
    #determine q value threshold or p value threshold
    q_use = TRUE
    while( sum(q_values<q_thresh) <= sig_gene )
    {
      q_thresh = q_thresh + 0.05
      if( q_thresh > 0.25)
      {
        q_use = FALSE
        break
      }
    }
    q_use = FALSE #use p value at last
    p_use = TRUE
    if( q_use == FALSE )
    {
      while( sum(p_values<p_thresh) <= sig_gene  )
      {
        p_thresh = p_thresh + 0.05
        if(p_thresh > 0.15)
        {
          p_use = FALSE
          break
        }
      }
    }
    
    #identify significant genes
    if( q_use == TRUE )
    {
      sig_ix = q_values<=q_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"q_thresh",q_thresh))
    }
    if( (q_use==F) && (p_use==T) )
    {
      sig_ix = p_values <= p_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"p_thresh",p_thresh))
    }
    if( (q_use==F) && (p_use==F) )
    {
      sig_ix = sort(p_values,decreasing=F,index.return=T)$ix
      sig_ix = sig_ix[1:sig_gene]
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"sig_gene",sig_gene))
    }
  }
  
  
  
  
}

find_0 <- function(train_dat,test_dat)
{
  #all row as genes and column as samples
  delete_ix = c()
  for(i in 1:nrow(train_dat))
  {
    if(sum(train_dat[i,]==0)==ncol(train_dat))
    {
      delete_ix = c(delete_ix, i)
    }
  }
  
  if( length(delete_ix>0) )
  {
    train_dat = train_dat[-delete_ix,]
    test_dat = test_dat[-delete_ix,]
    return(list(train_dat,test_dat))
  }
  if( length(delete_ix)==0 )
  {
    return(list(train_dat,test_dat))
  }
  
}

cancer_dummy <- function( data,info )
{
  #data: genes as rows and samples as columns
  
  #identify cancer types of each patient
  cancer_type = as.character( info$cancer[match(colnames(data),info$patient)] )
  cancer_type = as.factor(cancer_type)
  
  #turn cancer type to dummy variables
  tmp_df = data.frame(cancer_type,rep(1,times=length(cancer_type)))
  colnames(tmp_df) = c("cancer","tmp")
  cancers = t(model.matrix(tmp~cancer,tmp_df)[,-1])
  colnames(cancers) = colnames(data)
  row_names = substring(rownames(cancers),7)
  rownames(cancers) = row_names
  
  data = rbind(data,cancers)
  
  return(list(data,row_names))
  
}


dummy_to_test <- function(data, info, dummy)
{
  #data: row as features, col as patients
  #dummy: cancer type as dummy variables
  
  dummy_mat = matrix(0,nrow=length(dummy),ncol=ncol(data))
  pats = colnames(data)
  for(i in 1:ncol(data))
  {
    curr_pat = pats[i]
    curr_cancer = as.character( info$cancer[match(curr_pat,as.character(info$patient))] )
    ix = match(curr_cancer,as.character(dummy))
    
    if( !is.na(ix) )
    {
      dummy_mat[ix,i] = 1
    }
    
  }
  
  dimnames(dummy_mat) = list(as.character(dummy),colnames(data))
  data = rbind(data,dummy_mat)
  
  return(data)
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
#cluster
args <- commandArgs(trailingOnly=TRUE)
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")
test_fold = as.numeric(as.character(args[3]))#the current test fold, ranging from 1 to 5
output_folder = args[4] #no "/" at the end
freq = as.numeric(args[5]) #should be a adjustble parameter

#desktop
# setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv")
# cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
# cisplatin.info = read.table("cisplatin.gistic2.5_fold_cv.mat.txt",sep="\t",header=T,quote="")
# test_fold=1
# output_folder = "."
# freq = 0.8 
#both
test_fold = test_fold + 3
test_type = "regression"
#for debug
# test_type = "wilcox"


#libraries
library(glmnet)
library(doParallel)
library(foreach)
library(pracma)
no_cores = detectCores()

###preprocess data###
core.info = cisplatin.info
core.cancer = c("CESC","LUAD", "BLCA")

##find data##
#test data
test.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="validation"])
test.info = core.info[as.character(core.info[,test_fold])=="validation",]
test.dat = as.matrix(cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))])
test.resp = as.character(test.info$response[match(test.pats,as.character(test.info$patient))])
#for debug
#test.resp = as.character(test.info$cancer[match(test.pats,as.character(test.info$patient))])
#test.resp[test.resp!="CESC"] = "Others"
#

#train data
train.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="train"])
train.pats = sample(train.pats,size=length(train.pats),replace=FALSE)
train.info = core.info[as.character(core.info[,test_fold])=="train",]
train.dat = as.matrix(cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))])

#delete all zero genes
# tmp_list = find_0(train.dat,test.dat)
# train.dat = tmp_list[[1]]
# test.dat = tmp_list[[2]]

##select differential genes across cancer types##
# cl = makeCluster(no_cores)
# registerDoParallel(cl)
# list_tmp = test_gene(train.dat, test.dat, cisplatin.info, 
#                      type=test_type, parallel=T)
# stopImplicitCluster()
# stopCluster(cl)
# train.dat = list_tmp[[1]]
# test.dat = list_tmp[[2]]
# p_type = list_tmp[[3]]
# thresh = list_tmp[[4]]
# gene_num = nrow(train.dat)
# tmp_str = paste("With ",p_type," and threshold ",thresh,", ",
#                 gene_num," genes are remained",sep="")
# print(tmp_str)


#add cancer type#
gene_num = nrow(train.dat)
tmp_list = cancer_dummy( train.dat, cisplatin.info)
train.dat = tmp_list[[1]]
dummy = tmp_list[[2]]
test.dat = dummy_to_test(test.dat, cisplatin.info, dummy)


#delete all molecular data
train.dat = train.dat[-seq(1,gene_num,by=1),]
test.dat = test.dat[-seq(1,gene_num,by=1),]

###estimated robust features###
#basic parameters#
BS = 1000
alphas = seq(0.1,1,by=0.1)

#bootstrap samples#
#train patients response
train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
list.bs_mat = bootstrap_sample( colnames(train.dat), train.resp, BS=BS)
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]

##estimates recurrent features##
#model training#

best_lambda = matrix(NA,nrow=length(alphas),ncol=BS)
best_auc = matrix(NA,nrow=length(alphas),ncol=BS)
for(i in 1:length(alphas))
{
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  train_list <- foreach( bs=1:BS, .packages="glmnet",
                         .combine='comb', .multicombine=TRUE,
                         .init=list(list(), list()) ) %dopar%
  {
    #bootstrap sample data
    curr.train_dat = train.dat[ ,match(bs_mat.pats[,bs],colnames(train.dat)) ]
    curr.train_resp = bs_mat.resp[,bs]
  
    #record the best models in each bootstrap sample
    cv_fit = cv.glmnet( t(curr.train_dat), as.factor(curr.train_resp), 
                        family="binomial", type.measure="auc" )
  
    ix = match(cv_fit$lambda.1se,cv_fit$lambda)
    #beta = cv_fit$glmnet.fit$beta[,ix]
    auc = cv_fit$cvm[ix]
    return( list(cv_fit$lambda.1se,auc) )
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  curr_lambda = unlist(train_list[[1]])
  curr_auc = unlist( train_list[[2]] )
  
  best_lambda[i,] = curr_lambda
  best_auc[i,] = curr_auc
}

#find best models and corresponding features
cl = makeCluster(no_cores)
registerDoParallel(cl)
best_beta <- foreach(bs=1:BS,.combine='cbind',
                     .packages="glmnet") %dopar%
{
  #best model
  ix = which.max(best_auc[,bs])
  alpha = alphas[ix]
  lambda = best_lambda[ix,bs]
  
  #bootstrap data
  curr.train_dat = train.dat[ ,match(bs_mat.pats[,bs],colnames(train.dat)) ]
  curr.train_resp = bs_mat.resp[,bs]
  
  #train model
  glm_fit = glmnet( t(curr.train_dat), as.factor(curr.train_resp), 
                    family="binomial",alpha=alpha,lambda=lambda)
  
  curr_beta = as.vector(glm_fit$beta)
  
  return(curr_beta)
}
stopImplicitCluster()
stopCluster(cl)


#find recurrent features#
rownames(best_beta) = rownames(train.dat)
list_features = gene_selection(t(best_beta),freq=freq)
while( length(list_features[[2]]) <= 1 )
{
  freq = freq - 0.05
  list_features = gene_selection(t(best_beta),freq=freq)
}
tmp_str = paste( "Frequency ",freq," to select recurrent features" )
print(tmp_str)

#refine data by recurrent features
train.dat = train.dat[list_features[[2]],]
test.dat = test.dat[list_features[[2]],]
feature_freq = list_features[[3]]


##refit models and test##
#bootstrap samples#
list.bs_mat = bootstrap_sample( colnames(train.dat), train.resp, BS=BS)
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]

#fit models and test
cl = makeCluster(no_cores)
registerDoParallel(cl)
test_list <- foreach(bs=1:BS,.packages="glmnet") %dopar%
{
  #bootstrap sample
  curr.train_dat = train.dat[ ,match(bs_mat.pats[,bs],colnames(train.dat)) ]
  curr.train_resp = as.factor(bs_mat.resp[,bs])
  
  #train best model by lasso logistic
  #cv_fit = cv.glmnet( t(curr.train_dat), as.factor(curr.train_resp), 
  #                    family="binomial",type.measure="auc" )
  #fit the best model,should use the logistic regression
  glm_fit = glmnet(x=t(curr.train_dat),y=as.factor(curr.train_resp),
               lambda=0, family="binomial")
  
  #fit logistic
  #glm_fit = glm.fit(x=t(curr.train_dat),y=as.factor(curr.train_resp), 
  #                  family=binomial)
  #   df = data.frame(t(curr.train_dat),curr.train_resp)
  #   colnames(df) = c(rownames(curr.train_dat),"response")
  #   glm_fit = glm(response~., data=df, family=binomial())

  #predict on test data
  pred = predict( object=glm_fit,newx=t(test.dat), type="response" )

  return(as.vector(pred))
  
}
stopImplicitCluster()
stopCluster(cl)

test_score = do.call(cbind,test_list)
rownames(test_score) = test.pats


###output results###
#feature frequncy#
tmp_str = paste(output_folder,"/pan.elanet.feature_freq.test_",test_fold-3,".20150701.tiff",sep="")
tiff(tmp_str)
hist(feature_freq,main="Frequency of hittring for genes",xlab="hitting Freq",ylab="Freq",50)
dev.off()

##final results##
#selected features#
tmp_str = paste(output_folder,"/pan.elanet.feature.test_",test_fold-3,".20150701.txt",sep="")
write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")


#test by each model#
tmp_str = paste(output_folder,"/pan.elanet.mid_res.test_",test_fold-3,".20150701.txt",sep="")
write.table(t(test_score),tmp_str,col.names=T,row.names=T,sep="\t",quote=F)


##plot TEST performance##
tmp_str = paste(output_folder,"/pan.elanat.test_",test_fold-3,".20150708.tiff",sep="")
tiff(tmp_str)
#roc
roc = ensemble_roc(t(test_score),test.resp,"sensitive")
plot(roc[,2],roc[,1],"b",cex=0.8,xlab="FPR",ylab="TPR",
     xlim=c(0,1),ylim=c(0,1))
auc = trapz(roc[,2],roc[,1])
tmp_str = paste("AUC = ",round(auc,digits=2),sep="")
title("Test performance on Pan-Caner",tmp_str)
dev.off()


####THE END####


