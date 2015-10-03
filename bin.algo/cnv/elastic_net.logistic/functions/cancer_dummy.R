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