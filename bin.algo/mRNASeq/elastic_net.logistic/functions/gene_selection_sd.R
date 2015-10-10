gene_selection_sd <- function(data,sd=1.5)
{
  #identify genomic signature
  #data: row as genes, column as bootstap 1~100
  
  #the "single side frequency"
  freq = data
  freq[freq>0] = 1
  freq[freq<0] = -1
  side_freq = rowSums(freq)
  side_freq = side_freq/ncol(data)
  
  #the average weights
  weights = rowMeans(data)
  
  mean.freq = mean(side_freq)
  mean.weight = mean(weights)
  sd.freq = sd(side_freq)
  sd.weight = sd(weights)
  
  #identify genes
  ix.freq = c(which( side_freq>=(mean.freq+sd*sd.freq) ),
              which(side_freq<=(mean.freq-sd*sd.freq)) )
  ix.weight = c(which( weights>=(mean.weight+sd*sd.weight) ),
                which( weights<=(mean.weight-sd*sd.weight)) )
  
  ix = intersect(ix.freq,ix.weight)
  genes = rownames(data)[ix]
  
  freq_weight = cbind(c(side_freq[ix],side_freq[-ix]),c(weights[ix],weights[-ix]))
  rownames(freq_weight) = c(genes,rownames(data)[-ix])
  return( list(genes,ix,freq_weight) )
 
}