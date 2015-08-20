####for mRNAseq/expression data###
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

filter_mRNA <- function( data, train_pats, test_pats, low.thresh="Q1", type="dispersion" )                        
{
  #input: gene as rows and samples as columns
  #output: same as input
  #data include both train and test
  #low.thresh is the iqr/median threshold to filte lowly expressed genes, only for "dispersion" type
  #if type == "IQR", then genes with IQR smaller than Q1 will be filtered out
  
  #get data
  train_dat.ix = match( train_pats,colnames(data) )
  train.dat = data[,train_dat.ix]
  test_dat.ix = match(test_pats,colnames(data))
  test.dat = data[,test_dat.ix]
  
  #filter genes with low IQR/median
  #in train data, calculate coefficient of dispersion
  train_dat = as.matrix( t(train.dat) )
  if( type == "IQR")
  {
    keep_ix = vector(length=ncol(train_dat),mode="logical")
    for(i in 1:ncol(train_dat))
    {
      
      iqr = IQR(train_dat[,i],na.rm=T)
      q1 = quantile(train_dat[,i],na.rm=T)[2]
      if( iqr>= q1)
      {
        keep_ix[i] = TRUE
      }
      if( iqr < q1)
      {
        keep_ix[i] = FALSE
      }
    }
    
  }
  if( type == "dispersion")
  {
    #calculate coefficient of dispersion
    iqr = vector(length=ncol(train_dat),mode="numeric")
    for(i in 1:ncol(train_dat))
    {
      curr_quantile = quantile(train_dat[,i],na.rm=T)
      tmp_1 = curr_quantile[4] - curr_quantile[2]
      tmp_2 = curr_quantile[4] + curr_quantile[2]
      iqr[i] = tmp_1/tmp_2
    }
    #identify threshold and genes to be deleted
    thresh = quantile(iqr)[2]
    keep_ix = iqr>=thresh
  }
  
  
  
  #filter data
  train.dat = as.matrix(train.dat[keep_ix,])
  test.dat = as.matrix(test.dat[keep_ix,])
  
  return( list(train.dat,test.dat) )
  
}


exp_norm <- function(data, train_pats, test_pats, info)
{
  ##input
  #data with each row as genes, each column as samples, rownames
  ##output
  #data with each row as genes, each column as samples;output is matrice
  #NOTICE:  test dat normalization use the mean and sd of train data, RIGHT or WRONG?
  ##method
  #normalize gene expression of each cancer type
  #the default normalization method is z-score
  
  #find data
  train_ix = match(train_pats, colnames(data) )
  test_ix = match(test_pats, colnames(data) )
  train.dat = as.matrix(data[,train_ix])
  test.dat = as.matrix(data[,test_ix])
  
  #get cancer
  info_ix = match( unique(train_pats),as.character(info$patient) )
  train_cancer = unique(as.character(info$cancer[info_ix]))
  info_ix = match( unique(test_pats),as.character(info$patient) )
  test_cancer = unique(as.character(info$cancer[info_ix]))
  all_cancer = union(train_cancer,test_cancer)
  
  #train sample number in each cancer
  info_ix = match( unique(train_pats),as.character(info$patient) )
  train_cancer_table = table(as.character(info$cancer[info_ix]))
  cancer_num = vector(mode="numeric",length=length(all_cancer))
  for(i in 1:length(train_cancer_table))
  {
    ix = match(names(train_cancer_table)[i],all_cancer)
    cancer_num[ix] = train_cancer_table[i]
  }
  
  #mean values and sd values of each gene in each cancer
  mean_nums = matrix(NA,nrow=nrow(data),ncol=length(all_cancer))
  dimnames( mean_nums ) = list(rownames(data),all_cancer)
  sd_nums = matrix(NA,nrow=nrow(data),ncol=length(all_cancer))
  dimnames(sd_nums) = list(rownames(data),all_cancer)
  for(ca in 1:length(all_cancer))
  {
    #if in one cancer, no more than 1 sample
    if(cancer_num[ca]<=1)
    {
      mean_nums[,ca] = NA
      sd_nums[,ca] = NA
    }
    if( cancer_num[ca]> 1)
    {
      #find patients of train data in this cancer
      curr_pats = as.character( info$patient[as.character(info$cancer)==all_cancer[ca]] )
      curr_ix = c()
      dat_pats = colnames(train.dat)
      for( i in 1:ncol(train.dat) )
      {
        if(dat_pats[i] %in% curr_pats  )
        {
          curr_ix = c(curr_ix,i)
        }
      }
      curr_dat = train.dat[,curr_ix]
      
      #normalize train data of current cancer
      for(i in 1:nrow(data))
      {
        mean_nums[i,ca] = mean( curr_dat[i,], na.rm=T )
        sd_nums[i,ca] = sd( curr_dat[i,],na.rm=T )
        train.dat[i,curr_ix] = ( train.dat[i,curr_ix] - mean_nums[i,ca]  ) / sd_nums[i,ca]
      }
      
      #normalize test data of current cancer
      curr_test_ix = c()
      dat_pats = colnames(test.dat)
      for( i in 1:ncol(test.dat)  )
      {
        if( dat_pats[i] %in% curr_pats)
        {
          curr_test_ix = c(curr_test_ix,i)
        }
      }
      if(length(curr_test_ix)>0)
      {
        for(i in 1:nrow(test.dat))
        {
          test.dat[i,curr_test_ix] = ( test.dat[i,curr_test_ix] - mean_nums[i,ca] ) / sd_nums[i,ca]
        }
      }
    }
    
  }
  
  #normalization cancer existed in test data
  #but not in train data
  na_cancer.ix = which(cancer_num==0)
  if( length(na_cancer.ix)>0 )
  {
    #estimated mean value and sd value of na cancers by other cancer
    mean_nums[,na_cancer.ix] = rowMeans(mean_nums,na.rm=T)
    mean_na = rowMeans(mean_nums,na.rm=T)
    sd_nums[,na_cancer.ix] = rowMeans(sd_nums,na.rm=T)
    sd_na = rowMeans(sd_nums,na.rm=T)
    
    #missing cancers
    na_cancers = all_cancer[na_cancer.ix]
    
    #find patients in test data
    dat_pats = colnames(test.dat)
    for(i in 1:length(dat_pats))
    {
      curr_cancer = as.character(info$cancer[match(dat_pats[i],as.character(info$patient))])
      if( curr_cancer %in% na_cancers )
      {
        test.dat[,i] = (test.dat[,i] - mean_na) / sd_na
      }
    }
    
  }
  
  return(list(train.dat,test.dat))
  
}


test_gene <- function(train_dat, test_dat, info, sig_gene = 500,
                      p_thresh=0.05,q_thresh=0.05, type="ttest")
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
    for( i in 1:nrow(train_dat) )
    {
      sen_dat = as.numeric(train_dat[i,sen_ix])
      insen_dat = as.numeric(train_dat[i,insen_ix])
      test_fit = t.test(sen_dat,insen_dat,na.action="na.omit")
      p_values[i] = test_fit$p.value
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
    for( i in 1:nrow(train_dat) )
    {
      sen_dat = as.numeric(as.character(train_dat[i,sen_ix]))
      insen_dat = as.numeric(as.character(train_dat[i,insen_ix]))
      test_fit = wilcox.test(sen_dat,insen_dat,na.action="na.omit")
      p_values[i] = test_fit$p.value
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
    p_values = vector( length=nrow(train_dat),mode="numeric" )
    
    #perform test
    for(i in 1:nrow(train_dat) )
    {
      fit = glm(responses~as.numeric(as.character(train_dat[i,]))+cancers,family=binomial)
      sum_fit = summary(fit)
      coeff = sum_fit$coefficients
      p_values[i] = coeff[2,ncol(coeff)]
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


impute_NA <- function(train_dat,test_dat,gene_num)
{
  #input data: row as genes and col as samples
  #output data:row as genes and col as samples
  #usage: before glmnet function
  curr_train = train_dat[1:gene_num,]
  curr_test = test_dat[1:gene_num,]
  
  mean_nums = rowMeans( curr_train,na.rm=T )
  for(i in 1:nrow(train_dat))
  {
    train_na = which( is.na(train_dat[i,]) )
    test_na = which(is.na(test_dat[i,]))
    
    if( length(train_na)>0)
    {
      curr_train[i,train_na] = mean_nums[i]
    }
    if( length(test_na)>0)
    {
      curr_test[i,test_na] = mean_nums[i]
    }
  }
  
  train_dat[1:gene_num,] = curr_train
  test_dat[1:gene_num,] = curr_test
  
  return(list(train_dat,test_dat))
  
}


gene_selection <- function(data,freq=0.8)
{
  #data: row as BS and column as genes
  new_data = data; new_data[new_data!=0] = 1
  
  feature_freq = colSums(new_data)
  
  feature_freq = feature_freq/nrow(data)
  
  rownames(feature_freq) = rownames(data)
  
  ix = which(feature_freq>=freq)
  
  genes = colnames(data)[ix]
  
  return(list(genes,ix,feature_freq))
  
}


ensemble_roc <- function(score_mat,class,target_class,step=0.0001)
{
  #score_mat: row as one model, column as patients
  #class: responses
  #return value: tpr and fpr under each cut
  #cut: seq(0,1,by=0.01)
  
  class_1 = target_class
  class_2 = setdiff( unique(class), target_class )
  if( sum(class==class_1) > sum(class==class_2) )
  {
    class_p = class_2
    class_n = class_1
  }
  if( sum(class==class_1) <= sum(class==class_2) )
  {
    class_p = class_1
    class_n = class_2
  }
  
  
  
  cut_vec = seq(0,1,by=step)
  tpr_vec = vector(length=length(cut),mode="numeric")
  fpr_vec = vector(length=length(cut),mode="numeric")
  
  for( i in 1:length(cut_vec) )
  {
    curr_mat = matrix( NA,nrow=nrow(score_mat),ncol=ncol(score_mat) )
    curr_mat[ score_mat>=cut_vec[i] ] = 1
    curr_mat[ score_mat< cut_vec[i] ] = 0
    
    curr_score = colMeans( curr_mat )
    
    curr_class = vector(length=length(curr_score),mode="character")
    curr_class[curr_score>=0.5] = class_1
    curr_class[curr_score<0.5] = class_2
    
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


ensemble_pr <- function(score_mat,class,target_class,step=0.0001)
{
  #score_mat: row as one model, column as patients
  #class: responses
  #return value: tpr and fpr under each cut
  #cut: seq(0,1,by=0.01)
  
  class_1 = target_class
  class_2 = setdiff( unique(class), target_class )
  if( sum(class==class_1) > sum(class==class_2) )
  {
    class_p = class_2
    class_n = class_1
  }
  if( sum(class==class_1) <= sum(class==class_2) )
  {
    class_p = class_1
    class_n = class_2
  }
  
  
  
  cut_vec = seq(0,1,by=step)
  prec_vec = vector(length=length(cut),mode="numeric")
  rec_vec = vector(length=length(cut),mode="numeric")
  
  
  
  
  for( i in 1:length(cut_vec) )
  {
    curr_mat = matrix( NA,nrow=nrow(score_mat),ncol=ncol(score_mat) )
    curr_mat[ score_mat>=cut_vec[i] ] = 1
    curr_mat[ score_mat< cut_vec[i] ] = 0
    
    curr_score = colMeans( curr_mat )
    
    curr_class = vector(length=length(curr_score),mode="character")
    curr_class[curr_score>=0.5] = class_1
    curr_class[curr_score<0.5] = class_2
    
    #mytable = table( curr_class,class )
    mytable = matrix(0,nrow=2,ncol=2)
    dimnames(mytable) = list( c(class_p,class_n), c(class_p,class_n) )
    
    mytable[1,1] = sum(curr_class[class==class_p]==class_p)
    mytable[2,1] = sum(curr_class[class==class_p]==class_n)
    mytable[1,2] = sum(curr_class[class==class_n]==class_p)
    mytable[2,2] = sum(curr_class[class==class_n]==class_n)
    
    if( sum(mytable[,1])!=0 )
    {
      rec_vec[i] = mytable[1,1]/sum(mytable[,1])
    }
    if( sum(mytable[,1])==0 )
    {
      rec_vec[i] = NA 
    }
    if( sum(mytable[1,])!=0 )
    {
      prec_vec[i] = mytable[1,1]/sum(mytable[1,])
    }
    if( sum(mytable[1,])==0 )
    {
      tpr_vec[i] = NA
    }
    
  }
  
  pr_mat = cbind(prec_vec,rec_vec)
  colnames(pr_mat) = c("precision","recall")
  return(pr_mat)
}


####load data####
#cluster
args <- commandArgs(trailingOnly=TRUE)
cisplatin.dat = read.table(args[1],header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(args[2],sep="\t",header=T,quote="")
test_fold = as.numeric(as.character(args[3]))
core.cancer = as.character(args[4])
output_folder = args[5] #no "/" at the end
#desktop
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/mRNAseq")
cisplatin.dat = read.table("cisplatin.mRNAseq.gdac_20141206.preprocess.txt",header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table("cisplatin.mRNAseq_fold_cv.mat.txt",sep="\t",header=T,quote="")
core.cancer = "CESC"
test_fold=1
output_folder = "."
#both
test_fold = test_fold + 3

#libraries
library(glmnet)
library(PRROC)
library(ROCR)
library(doParallel)
library(foreach)
library(pracma)
no_cores = detectCores()

###preprocess data###
##find data##
core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer,]
#test data
test.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="validation"])
test.info = core.info[as.character(core.info[,test_fold])=="validation",]

#test patients response
test.resp = as.character(test.info$response[match(test.pats,as.character(test.info$patient))])
test.resp.lab = vector(mode="numeric",length=length(test.resp))
test.resp.lab[which(test.resp=="insensitive")] = 1
test.resp.lab[which(test.resp=="sensitive")] = 0

#train data
train.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="train"])
train.pats = sample(train.pats,size=length(train.pats),replace=FALSE)
train.info = core.info[as.character(core.info[,test_fold])=="train",]


##filter lowly expressed genes##
data.tmp = filter_mRNA(cisplatin.dat, train.pats, test.pats, 
                       low.thresh="Q1", type="dispersion" )
train.dat = data.tmp[[1]]; test.dat = data.tmp[[2]]


##normalization in each cancer##
data.tmp = exp_norm(cbind(train.dat,test.dat), train.pats, test.pats, cisplatin.info)
train.dat = data.tmp[[1]]
test.dat = data.tmp[[2]]

##select differential genes across cancer types##
list_tmp = test_gene(train.dat, test.dat, cisplatin.info, type="ttest")
train.dat = list_tmp[[1]]
test.dat = list_tmp[[2]]
test_type = list_tmp[[3]]
thresh = list_tmp[[4]]
tmp_str = paste("With ",test_type," and threshold ",thresh,", ",
                nrow(train.dat)," genes are remained",sep="")
print(tmp_str)

#impute NAs#
tmp_list = impute_NA(train.dat,test.dat,nrow(train.dat))
train.dat = tmp_list[[1]]
test.dat = tmp_list[[2]]

###estimated robust features###
#basic parameters#
BS = 100
folds = 5
alphas = seq(0.1,1,by=0.1)
#bootstrap samples#
#train patients response
train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
train.resp.lab = vector(mode="numeric",length=length(train.resp))
train.resp.lab[which(train.resp=="insensitive")] = 1
train.resp.lab[which(train.resp=="sensitive")] = 0

list.bs_mat = bootstrap_sample( colnames(train.dat), train.resp, BS=BS)
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]


##estimates recurrent features##
#model training#
best_params = matrix(NA,nrow=BS,ncol=4)
colnames(best_params) = c("alpha","lambda","AUC","nzero")
best_betas = matrix(NA,nrow=BS,ncol=nrow(train.dat))
colnames(best_betas) = rownames(train.dat)
for( bs in 1:BS  )
{
  curr.train_dat = train.dat[ ,match(bs_mat.pats[,bs],colnames(train.dat)) ]
  curr.train_resp = bs_mat.resp[,bs]
  
  
  #train models for current bootstrap sample
  bs_params = matrix( NA,nrow=length(alphas),ncol=4 )
  colnames(bs_params) = c("alpha","lambda","AUC","nzero")
  bs_betas = matrix( NA, nrow=length(alphas), ncol=nrow(curr.train_dat) )
  colnames(bs_betas) = rownames(curr.train_dat)
  for( i in 1:length(alphas) )
  {
    #do parallel cross validation
    cl = makeCluster(no_cores)
    registerDoParallel(cl)
    cv_fit = cv.glmnet( t(curr.train_dat), as.factor(curr.train_resp), family="binomial",
                        type.measure="auc",alpha=alphas[i], parallel=TRUE,nfolds=3)
    stopImplicitCluster()
    stopCluster(cl)
    
    #record current best model,under specific alpha
    bs_params[i,1] = alphas[i]
    bs_params[i,2] = cv_fit$lambda.1se
    bs_params[i,3] = cv_fit$cvm[match(cv_fit$lambda.1se,cv_fit$lambda)]
    bs_params[i,4] = cv_fit$nzero[match(cv_fit$lambda.1se,cv_fit$lambda)]
    full_models = cv_fit$glmnet.fit 
    beta = full_models$beta[,match(cv_fit$lambda.1se,full_models$lambda)]
    bs_betas[i,] = beta[match(colnames(bs_betas),names(beta))]
  }
  
  #find best model for current bootstrap sample
  
  best_ix = which.max(bs_params[,3])
  best_params[bs,3] = bs_params[best_ix,3]
  best_params[bs,1] = bs_params[best_ix,1]
  best_params[bs,2] = bs_params[best_ix,2]
  best_params[bs,4] = bs_params[best_ix,4]
  
  best_betas[bs,] = bs_betas[best_ix,]
  
}

#recurrent features#
freq = 0.8
list_features = gene_selection(best_betas,freq=freq)
#why 10? may be incorrect
while( length(list_features[[2]]) <= 10 )
{
  freq = freq - 0.05
  list_features = gene_selection(best_betas,freq=freq)
}

#refine data by recurrent features
train.dat = train.dat[list_features[[2]],]
test.dat = test.dat[list_features[[2]],]
feature_freq = list_features[[3]]


##refit models and test##
#bootstrap samples#
list.bs_mat = bootstrap_sample( colnames(train.dat), train.resp, BS=BS)
bs_mat.pats = list.bs_mat[[1]]
bs_mat.resp = list.bs_mat[[2]]


cv_score = matrix( NA,nrow=BS,ncol=nrow(bs_mat.pats) )
cv_auc = vector(length=BS,mode="numeric")
#colnames(cv_score) = colnames(train.dat)
test_score = matrix( NA,nrow=BS,ncol=ncol(test.dat) )
colnames(test_score) = colnames(test.dat)
test_class = matrix( "NULL",nrow=BS,ncol=ncol(test.dat) )
colnames(test_class) = colnames(test.dat)
for(bs in 1:BS)
{
  #data
  curr.train_dat = train.dat[,match(bs_mat.pats[,bs],colnames(train.dat)) ]
  curr.train_resp = as.factor(bs_mat.resp[,bs])
  #find best model
  bs_params = matrix( NA,nrow=length(alphas),ncol=4 )
  colnames(bs_params) = c("alpha","lambda","AUC","nzero")
  bs_betas = matrix( NA, nrow=length(alphas), ncol=nrow(curr.train_dat) )
  colnames(bs_betas) = rownames(curr.train_dat)
  for( i in 1:length(alphas) )
  {
    #do parallel cross validation
    cl = makeCluster(no_cores)
    registerDoParallel(cl)
    cv_fit = cv.glmnet( t(curr.train_dat), as.factor(curr.train_resp), family="binomial",
                        type.measure="auc",alpha=alphas[i], parallel=TRUE,nfolds=3)
    stopImplicitCluster()
    stopCluster(cl)
    
    #record current best model,under specific alpha
    bs_params[i,1] = alphas[i]
    bs_params[i,2] = cv_fit$lambda.1se
    bs_params[i,3] = cv_fit$cvm[match(cv_fit$lambda.1se,cv_fit$lambda)]
    bs_params[i,4] = cv_fit$nzero[match(cv_fit$lambda.1se,cv_fit$lambda)]
    full_models = cv_fit$glmnet.fit 
    beta = full_models$beta[,match(cv_fit$lambda.1se,full_models$lambda)]
    bs_betas[i,] = beta[match(colnames(bs_betas),names(beta))]
  }
  
  best_ix = which.max(bs_params[,3])
  alpha_best = bs_params[best_ix,1]
  lambda_best = bs_params[best_ix,2]
  cv_auc[bs] = bs_params[best_ix,3]
  
  
  #fit the best model
  fit = glmnet(x=t(curr.train_dat),y=curr.train_resp,alpha=alpha_best,
               lambda=lambda_best, family="binomial")
  #this is wrong, even the prediction results 
  pred_cv = predict(object=fit,newx=t(curr.train_dat), family="response")
  
  cv_score[bs,] = as.vector(pred_cv)
  ##############################################################
  #notice: the cv performance needs a real cross-validation
  # the cutoff should be performed based on this cross-vlaidation
  #therefore,the CV performance is replaced by train error
  #therefore,the classification cutoff is 0.5
  ##############################################################
  
  #predict on test data
  pred = predict( object=fit,newx=t(test.dat), type="response" )
  test_score[bs,] = pred[ match( colnames(test_score), rownames(pred) ) ]
  
  #record predicted classes
  pred = predict( object=fit,newx=t(test.dat),type="class")
  test_class[bs,] = pred[ match( colnames(test_class), rownames(pred) ) ]
  
  
}


##plot TEST performance##
#draw ROC and precision-recall curves
tmp_str = paste(output_folder,"/single.elanat.test_performance.test_",test_fold-3,".20150701.pdf",sep="")
pdf(tmp_str)
#roc
roc = ensemble_roc(test_score,test.resp,"sensitive")
plot(roc[,2],roc[,1],"b",cex=0.8,xlab="FPR",ylab="TPR",
     xlim=c(0,1),ylim=c(0,1))
auc = trapz(roc[,2],roc[,1])
tmp_str = paste("AUC = ",auc,sep="")
title("Test performance on CESC",tmp_str)
#pr
pr= ensemble_pr(test_score,test.resp,"sensitive")
plot(pr[,2],pr[,1],"b",cex=0.8,xlab="recall",ylab="precision",
     xlim=c(0,1),ylim=c(0,1))
pr_tmp = pr[which(!is.na(pr[,1])),]
pr_tmp = pr_tmp[which(!is.na(pr_tmp[,2])),]
auc = trapz(pr_tmp[,2],pr_tmp[,1])
tmp_str = paste("AUC = ",auc,sep="")
title("Test performance on CESC",tmp_str)

dev.off()

###output results###
##middle results##
#CV AUCs#
tmp_str = paste(output_folder,"/single.elasnet.cv_auc.test_",test_fold-3,"20150701.tiff",sep="")
tiff(tmp_str)
boxplot(cv_auc,main="AUCs of each models of bagging",ylab="AUC")
dev.off()

#eature frequncy#
tmp_str = paste(output_folder,"/single.elanet.feature_freq.test_",test_fold-3,".20150701.tiff",sep="")
tiff(tmp_str)
hist(feature_freq,main="Frequency of hittring for genes",xlab="hitting Freq",ylab="Freq",50)
dev.off()

##test results##
#test by each model#
tmp_str = paste(output_folder,"/single.elanet.mid_res.test_",test_fold-3,".20150701.txt",sep="")
write.table(test_score,tmp_str,col.names=T,row.names=F,sep="\t",quote=F)


##selected features##
tmp_str = paste(output_folder,"/single.elanet.feature.test_",test_fold-3,".20150701.txt",sep="")
write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")



####THE END####
