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