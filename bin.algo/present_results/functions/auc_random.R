auc_random <- function( roc, test="wilcox" )
{
  #roc is a list, each element is a 2-column matrix
  #1st column as y and 2nd column as x
  library(pracma)
  
  auc = c()
  for(i in 1:length(roc))
  {
    curr_auc = trapz(roc[[i]][,2],roc[[i]][,1])
    auc = c(auc,curr_auc)
  }
  
  
  if(test=="wilcox")
  {
    wilcox = wilcox.test(auc, alternative="greater", mu=0.5)
    hist(auc,20,xlab="AUC",ylab="Frequency",main=NULL)
    title("Histogram of AUC",paste("Wilcox Rank Sum Test pvalue =",wilcox$p.value))
    return(wilcox)
  }
  if(test=="ttest" )
  {
    t_test = t.test(auc, alternative="greater", mu=0.5)
    hist(auc,20,xlab="AUC",ylab="Frequency",main=NULL)
    title("Histogram of AUC",paste("T Test pvalue =",t_test$p.value))
    return(t_test)
  }
  if(test=="ztest")
  {
    library(BSDA)
    z_test = z.test(auc,alternative="greater",mu=0.5,sigma.x=sd(auc))
    hist(auc,20,xlab="AUC",ylab="Frequency",main=NULL)
    title("Histogram of AUC",paste("Z Test pvalue =",z_test$p.value))
    return(z_test)
  }
  
  
}