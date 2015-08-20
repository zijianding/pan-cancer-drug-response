#functions
ensemble_roc <- function(score_mat,class,target_class,step=0.0001)
{
  #score_mat: row as one model, column as patients
  #class: responses
  #return value: tpr and fpr under each cut
  #cut: seq(0,1,by=0.01)
  #currently target class is treated as negative classs
  #which is wrong, later we will corrected this
  # for debug
  #     score_mat = t(test_score)
  #     class = test.resp
  #     target_class = "sensitive"
  #     step=0.0001
  
  #
  
  
  #here is the problem
  class_1 = target_class # large score
  class_2 = setdiff( unique(class), target_class )
  #   if( sum(class==class_1) > sum(class==class_2) )
  #   {
  #     class_p = class_2
  #     class_n = class_1
  #   }
  #   if( sum(class==class_1) <= sum(class==class_2) )
  #   {
  #     class_p = class_1
  #     class_n = class_2
  #   }
  class_p = class_2 #small score
  class_n = class_1 #large score
  
  
  
  cut_vec = seq(0,1+step,by=step)
  cut_vec[length(cut_vec)] = Inf
  tpr_vec = vector(length=length(cut_vec),mode="numeric")
  fpr_vec = vector(length=length(cut_vec),mode="numeric")
  
  for( i in 1:length(cut_vec) )
  {
    curr_mat = matrix( NA,nrow=nrow(score_mat),ncol=ncol(score_mat) )
    curr_mat[ score_mat>=cut_vec[i] ] = 1 #large score
    curr_mat[ score_mat< cut_vec[i] ] = 0 #small score
    
    curr_score = colMeans( curr_mat )
    
    curr_class = vector(length=length(curr_score),mode="character")
    curr_class[curr_score>0.5] = class_1
    curr_class[curr_score<=0.5] = class_2
    
    #mytable = table( curr_class,class )
    mytable = matrix(0,nrow=2,ncol=2)
    dimnames(mytable) = list( c(class_p,class_n), c(class_p,class_n) )
    
    mytable[1,1] = sum(curr_class[class==class_p]==class_p)
    mytable[2,1] = sum(curr_class[class==class_p]==class_n)
    mytable[1,2] = sum(curr_class[class==class_n]==class_p)
    mytable[2,2] = sum(curr_class[class==class_n]==class_n)
    
    
    if( sum(mytable[,1])!=0 )
    {
      tpr_vec[i] = mytable[1,1]/sum(mytable[,1])
    }
    if( sum(mytable[,1])==0 )
    {
      tpr_vec[i] = NA 
    }
    if( sum(mytable[,2])!=0 )
    {
      fpr_vec[i] = mytable[1,2]/sum(mytable[,2])
    }
    if( sum(mytable[,2])==0 )
    {
      fpr_vec[i] = NA
    }
    
  }
  
  roc_mat = cbind(tpr_vec,fpr_vec)
  colnames(roc_mat) = c("tpr","fpr")
  return(roc_mat)
  
}

#10 features#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/artificial_dat/")
#balanced data
roc10_balance = read.table("roc.501_pos.10_features.4_selected.txt",
                           header=F,sep="\t",quote="")
roc10_balance = as.matrix(roc10_balance  )


#imbalanced data
roc10_im = read.table("roc.200_pos.10_features.2_selected.txt",
                      header=F,sep="\t",quote="")

roc10_im = as.matrix(roc10_im)

#1000 features#
#balanced data
roc1000_bal = read.table("roc.480_pos.1000_feature.52_selected.txt",
                         header=F,sep="\t",quote="")
roc1000_bal = as.matrix(roc1000_bal)


#imbalanced data
roc1000_im = read.table("roc.209_pos.1000_feature.41_selected.txt",
                        header=F,sep="\t",quote="")
roc1000_im = as.matrix(roc1000_im)









plot(roc10_balance[,2],roc10_balance[,1],cex=0.8,xlab="FPR",ylab="TPR",
     xlim=c(0,1),ylim=c(0,1),col="red",lty=1,type="l")
lines( roc10_im[,2],roc10_balance[,1],cex=0.8,
       col="red",lty=2,type="l"  )
lines(  roc1000_bal[,2],roc1000_bal[,1],cex=0.8,
        col="blue",lty=1,type="l")
lines(roc1000_im[,2],roc1000_im[,1],cex=0.8,
      col="blue",lty=2,type="l")
lines(seq(0,1,0.01),seq(0,1,0.01),lty=2,col="gray")
legend("bottomright",col=c("red","red","blue","blue"),
       lty=c(1,2,1,2),
       legend=c("balance data;10 feautures","imbalance data;10 features","balance data;1000 features","imbalance data;1000 features"))
title("Linear Model on Artificial data","1000 samples:balance data has 500 positive samples;imbalance 800")



