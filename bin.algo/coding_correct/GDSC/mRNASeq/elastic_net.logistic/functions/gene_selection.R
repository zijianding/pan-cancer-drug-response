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