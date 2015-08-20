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