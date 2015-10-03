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