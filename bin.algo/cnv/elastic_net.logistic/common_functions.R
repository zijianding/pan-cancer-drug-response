###pan-cancer###
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


###mRNA data###
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
    #keep_ix = vector(length=ncol(train_dat),mode="logical")
    iqr = vector(length=ncol(train_dat),mode="numeric")
    for(i in 1:ncol(train_dat))
    {
      
      curr_iqr = IQR(train_dat[,i],na.rm=T)
      #       q1 = quantile(train_dat[,i],na.rm=T)[2]
      #       if( iqr>= q1)
      #       {
      #         keep_ix[i] = TRUE
      #       }
      #       if( iqr < q1)
      #       {
      #         keep_ix[i] = FALSE
      #       }
      iqr[i] = curr_iqr
    }
    thresh = quantile(iqr)[2]
    keep_ix = iqr>=thresh
    
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


###CNV data###
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


###data management###
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


###parallel###
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}



####feature pre-selection###
#for debug
# train_dat=train.dat; test_dat = test.dat;info=cisplatin.info 
# type="regression"; parallel=T
#
# 统一使用p value
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

#下面的test_gene是用来确认代码准确性的,但还是利用肿瘤数据
test_gene <- function(train_dat, test_dat, info, sig_gene = 50, core.cancer,
                      p_thresh=0.05,q_thresh=0.05, type="ttest",parallel=F)
{
  #input data: row as genes and col as samples
  #output data:row as genes and col as samples
  #output significant threshold: q_thresh and p_thresh
  #type: default is t test; 
  #wilcox rank sum test, and logistic regression are also provided: "wilcox", "regression"
  #sig_gene, when no significant genes according to p/q value threshold
  #   train_dat=train.dat; test_dat = test.dat;info=cisplatin.info 
  #   type="wilcox"; parallel=T;core.cancer = core.cancer; sig_gene=50;
  #   p_thresh = 0.5; q_thresh = 0.5; 
  
  
  if( type=="ttest" )
  {
    #perform on train data
    train_pats = colnames(train_dat)
    curr_info = info[match(train_pats,as.character(info$patient)),]
    sen_ix = which(as.character(curr_info$cancer)==core.cancer)
    insen_ix = which(as.character(curr_info$cancer)!=core.cancer)
    
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
    sen_ix = which(as.character(curr_info$cancer)==core.cancer)
    insen_ix = which(as.character(curr_info$cancer)!=core.cancer)
    
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

#下面的test_gene是用来确认代码准确性的，利用人造数据
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


###find recurrent gene###
#consider positive/negative
#NOT consider scale, like mean and sd
#THE FREQUENCY is based on the difference value
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


###performance###
#consider for bagging logistic regression
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
