identify_class <- function(score_mat,class,target_class,cutoff=0.5)
{
  #score_mat: column as patients
  #target class: probability larger than 0.5
  mat = matrix( 0,nrow=nrow(score_mat),ncol=ncol(score_mat) )
  mat[score_mat>=cutoff] = 1
  vec = colMeans(mat)
  
  class_1 = target_class
  class_2 = setdiff(unique(class),class_1)
  
  pred = vector(length=length(vec),mode="character")
  pred[vec>=0.5] = class_1
  pred[vec< 0.5] = class_2
  
  indict = pred==class
  
  df = data.frame(class,pred,indict)
  colnames(df) = c("true","pred","indicator")
  rownames(df) = colnames(score_mat)
  
  mytable = table(pred,class)
  fdr = mytable[1,2]/sum(mytable[1,])
  
  #return(list(df,fdr))
  return(df)
  
}