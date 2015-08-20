bootstrap_sample <- function(obs,class,BS){
  #this function only generates BS samples
  #in each sample, two classes are balanced
  #the sample size equals to the number of obs(+1/-1 due to odd-even)
  #return a list
  
  sample.num = length( obs )
  
  bs_mat.sample = matrix( "NULL",ncol=BS,nrow=(2*sample.num) )
  bs_mat.resp = matrix( "NULL",ncol=BS,nrow=(2*sample.num) )
  
  for(bs in 1:BS)
  {
    #sample positive observations
    sample_obs = sample(obs,size=sample.num,replace=T)
    sample_obs_resp = class[match(sample_obs,obs)]
    
    bs_mat.sample[,bs] = sample_obs
    bs_mat.resp[,bs] = sample_obs_resp
  }
  
  bs_mat = list(bs_mat.sample,bs_mat.resp)
  
  return(bs_mat)
}