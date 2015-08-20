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