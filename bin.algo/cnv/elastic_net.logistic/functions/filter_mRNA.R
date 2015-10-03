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
