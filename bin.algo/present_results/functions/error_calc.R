error_calc <- function(data, info)
{
  #data as a list of matrice
  #each matrix: row as samples, column as true class, prediciton,etc
  #info: all information of patients by the drug
  
  #debug
  #data = error_res; info = cisplatin.info
  test_num = length(data)
  
  res_list = list()
  all_cancer = c("ALL") 
  
  for(K in 1:test_num)
  {
    patients = rownames(data[[K]])
    cancers = as.character( info$cancer[ match(patients,as.character(info$patient)) ] )
    cancer_type = unique( cancers )
    
    sen_vec = vector( mode="numeric",length=length(cancer_type)+1 )
    insen_vec = vector( mode="numeric", length=length(cancer_type)+1 )
    
    #each cancer
    for(i in 1:length(cancer_type))
    {
      ix = which( cancers==cancer_type[i] )
      curr_dat = data[[K]][ix,]
      
      #specifity
      sen = sum(curr_dat$indicator[curr_dat$true=="sensitive"]==T)
      sen_all = sum(curr_dat$true=="sensitive")
      sen_wrong = sen/sen_all
      
      #sensitivity
      insen = sum(curr_dat$indicator[curr_dat$true=="insensitive"]==T)
      insen_all = sum(curr_dat$true=="insensitive")
      insen_wrong = insen/insen_all
      
      sen_vec[i] = sen_wrong
      insen_vec[i] = insen_wrong
    }
    
    #all cancer
    #specificity
    sen = sum(data[[K]]$indicator[data[[K]]$true=="sensitive"]==T)
    sen_all = sum(data[[K]]$true=="sensitive")
    sen_vec[length(sen_vec)] = sen/sen_all
    #sensitivity
    insen = sum(data[[K]]$indicator[data[[K]]$true=="insensitive"]==T)
    insen_all = sum(data[[K]]$true=="insensitive")
    insen_vec[length(insen_vec)] = insen/insen_all
    
    
    
    curr_res = cbind(sen_vec,insen_vec)
    rownames(curr_res) = c(cancer_type,"ALL")
    colnames(curr_res) = c("sen","insen")
    
    res_list[[K]] = curr_res
    all_cancer  = union(all_cancer,cancer_type)
    
    
  }
  
  #pool all results
  all_cancer = unique(all_cancer)
  final_res = matrix(NA,nrow=2*length(all_cancer),ncol=test_num)
  rownames(final_res) = rep(all_cancer,each=2)
  resp_vec = rep(c("Specifity","Sensitivity"),times=length(all_cancer))
  
  for(K in 1:test_num)
  {
    curr_res = res_list[[K]]
    cancer_type = rownames(curr_res)
    
    for( i in 1:nrow(curr_res) )
    {
      cancer_ix = match(cancer_type[i],rownames(final_res))
      final_res[cancer_ix,K] = curr_res[i,1]
      final_res[cancer_ix+1,K] = curr_res[i,2]
    }
  }
  
  final = data.frame(rownames(final_res),resp_vec,final_res[,1])
  colnames(final) = c("cancer","Type","value")
  for(K in 2:test_num)
  {
    tmp = data.frame(rownames(final_res),resp_vec,final_res[,K])
    colnames(tmp) = c("cancer","Type","value")
    final = rbind( final,tmp )
  }
  
  #boxplot
  library(ggplot2)
  ggplot( aes(y=value,x=cancer),data=final ) + 
    geom_boxplot(aes(fill=Type)) + 
    geom_point(position=position_dodge(width=0.75), aes(fill=Type)) +
    ggtitle("Performance in Pan- and Single Cancer") +
    xlab("cancer types") +
    theme_set(theme_gray(base_size=25))
  
}