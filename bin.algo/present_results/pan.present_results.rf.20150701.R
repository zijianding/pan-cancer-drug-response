###library###
library(doParallel)
library(foreach)
no_cores = detectCores()
###functions###
read_table <- function(x)
{
  return( read.table(x,sep="\t",header=T,quote="",row.names=1) )
}

find_name <- function(x)
{
  arr = unlist(strsplit(x,split="\\/",perl=T))
  return(arr[length(arr)])
}

identify_class <- function(score_vec,class,target_class,cutoff=0.5)
{
  #score_mat: column as patients
  #target class: probability larger than 0.5
  
  #mat = matrix( 0,nrow=nrow(score_mat),ncol=ncol(score_mat) )
  #mat[score_mat>=cutoff] = 1
  #vec = colMeans(mat)
  vec = score_vec
  
  class_1 = target_class
  class_2 = setdiff(unique(class),class_1)
  
  pred = vector(length=length(vec),mode="character")
  pred[vec>=0.5] = class_1
  pred[vec< 0.5] = class_2
  
  indict = pred==class
  
  df = data.frame(class,pred,indict)
  colnames(df) = c("true","pred","indicator")
  rownames(df) = names(score_vec)
  
  mytable = table(pred,class)
  fdr = mytable[1,2]/sum(mytable[1,])
  
  #return(list(df,fdr))
  return(df)
  
}

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


ensemble_roc <- function(score_vec,class,target_class,step=0.0001)
{
  #score_mat: row as one model, column as patients
  #class: responses
  #return value: tpr and fpr under each cut
  #cut: seq(0,1,by=0.01)
  #currently target class is treated as negative classs
  #which is wrong, later we will corrected this
  # for debug
  score_vec = pan_score[[1]]
  class = resp
  target_class = "insensitive"
  step=0.0001
  
  #
  
  
  #here is the problem
  #in RF, large score="insensitive"
  class_1 = target_class # large score
  class_2 = setdiff( unique(class), target_class )
  class_n = class_2 #small score
  class_p = class_1 #large score
  
  
  
  cut_vec = seq(0,1+step,by=step)
  cut_vec[length(cut_vec)] = Inf
  tpr_vec = vector(length=length(cut_vec),mode="numeric")
  fpr_vec = vector(length=length(cut_vec),mode="numeric")
  
  for( i in 1:length(cut_vec) )
  {
    #curr_mat = matrix( NA,nrow=length(score_mat),ncol=ncol(score_mat) )
    curr_vec = vector("numeric",length(score_vec))
    curr_vec[ score_vec>=cut_vec[i] ] = 1 #large score
    curr_vec[ score_vec< cut_vec[i] ] = 0 #small score
    
    curr_score = curr_vec
    
    curr_class = vector(length=length(curr_score),mode="character")
    curr_class[curr_score>0.5] = class_p
    curr_class[curr_score<=0.5] = class_n
    
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
  
  roc_mat = cbind(rev(tpr_vec),rev(fpr_vec))
  colnames(roc_mat) = c("tpr","fpr")
  return(roc_mat)
  
}


ave_roc <- function(data, type="vertical",title_str,measure="quantile")
{
  #data as a list of matrice
  #each matrix: tpr/fpr
  #return: an average curve with error bar(sd)
  #for debug
  #       data = roc_res
  #       type= "vertical"
  #       measure = "quantile"
  #       title_str = "test"
  #now measure should be used along with quantila
  #threshold is currently not be used along with quantile
  
  
  library(pracma)
  data_num = length(data)
  
  #by vertical
  if( type=="vertical" )
  {
    #fpr points
    fpr_fix = seq(0,1,by=0.01)
    fpr_fix.sd_ix = seq(11,91,by=10)
    #tpr points
    tpr_fix = vector(length=length(fpr_fix),
                     mode="numeric")
    tpr_fix.sd = vector(mode="numeric",
                        length=length(fpr_fix.sd_ix))
    
    #calculate fpr
    fpr_fix.mat = matrix(NA,nrow=length(fpr_fix),
                         ncol=data_num)
    for( i in 1:data_num )
    {
      for( j in 2:(length(fpr_fix)-1) )
      {
        ix = which(data[[i]][,2]==fpr_fix[j])
        if( length(ix)>0 )
        {
          fpr_fix.mat[j,i] = max(data[[i]][ix,1])
        }
        if( length(ix)==0)
        {
          val = abs( data[[i]][,2]-fpr_fix[j] )
          ix = which( val == min(val) )
          curr_val = unique(data[[i]][ix,2])
          if( length(curr_val)>1)
          {
            x1 = data[[i]][max(which(data[[i]][,2]==curr_val[1])),2]
            y1 = data[[i]][max(which(data[[i]][,2]==curr_val[1])),1]
            x2 = data[[i]][min(which(data[[i]][,2]==curr_val[2])),2]
            y2 = data[[i]][min(which(data[[i]][,2]==curr_val[2])),1]
            y = (y2-y1)/(x2-x1)*fpr_fix[j]
            y = y + y2 - (y2-y1)/(x2-x1)*x2
            fpr_fix.mat[j,i] = y
            
          }
          if( length(curr_val)==1)
          {
            if(curr_val<fpr_fix[j])
            {
              if( length(ix) > 1)
              {
                ix = max(ix)
              }
              y1 = data[[i]][ix,1]
              x1 = data[[i]][ix,2]
              y2 = data[[i]][ix+1,1]
              x2 = data[[i]][ix+1,2]
              y = (y2-y1)/(x2-x1)*fpr_fix[j]
              y = y + y2 - (y2-y1)/(x2-x1)*x2
              fpr_fix.mat[j,i] = y
            }
            if( curr_val>fpr_fix[j] )
            {
              if( length(ix) > 1)
              {
                ix = min(ix)
              }
              y1 = data[[i]][ix-1,1]
              x1 = data[[i]][ix-1,2]
              y2 = data[[i]][ix,1]
              x2 = data[[i]][ix,2]
              y = (y2-y1)/(x2-x1)*fpr_fix[j]
              y = y + y2 - (y2-y1)/(x2-x1)*x2
              fpr_fix.mat[j,i] = y
            }
          }
          
          
          
        }
      }
      
    }
    fpr_fix.mat[1,] = 0
    fpr_fix.mat[nrow(fpr_fix.mat),] = 1
    
    tpr_fix = rowMeans( fpr_fix.mat )
    for( i in 1:length(fpr_fix.sd_ix) )
    {
      tpr_fix.sd[i] = sd(fpr_fix.mat[fpr_fix.sd_ix[i],])
    }
    
    tpr_quantile = t( apply(fpr_fix.mat,1,quantile) )
    
    
    #draw the average curve
    if( measure=="sd")
    {
      plot(fpr_fix, tpr_fix, "l",
           xlim=c(0,1), ylim=c(0,1),
           xlab="False Positive Rate",
           ylab="True Positive Rage",
           main=title_str)
      lines(seq(0,1,0.1),seq(0,1,0.1),lty=2,col="gray")
      segments(fpr_fix[fpr_fix.sd_ix],tpr_fix[fpr_fix.sd_ix]-tpr_fix.sd,
               fpr_fix[fpr_fix.sd_ix],tpr_fix[fpr_fix.sd_ix]+tpr_fix.sd)
      epsilon = 0.01
      segments(fpr_fix[fpr_fix.sd_ix]-epsilon,tpr_fix[fpr_fix.sd_ix]-tpr_fix.sd,
               fpr_fix[fpr_fix.sd_ix]+epsilon,tpr_fix[fpr_fix.sd_ix]-tpr_fix.sd)
      segments(fpr_fix[fpr_fix.sd_ix]-epsilon,tpr_fix[fpr_fix.sd_ix]+tpr_fix.sd,
               fpr_fix[fpr_fix.sd_ix]+epsilon,tpr_fix[fpr_fix.sd_ix]+tpr_fix.sd)
      return( list(cbind(tpr_fix,fpr_fix),cbind(fpr_fix.sd_ix,tpr_fix.sd)) )
    }
    
    if(measure=="quantile")
    {
      plot(fpr_fix,tpr_quantile[,3],"l",
           xlim=c(0,1),ylim=c(0,1),
           xlab="False Positive Rate",
           ylab="True Positive Rage",
           main=title_str)
      lines(seq(0,1,0.1),seq(0,1,0.1),lty=2,col="gray")
      
      segments(fpr_fix[fpr_fix.sd_ix],
               tpr_quantile[fpr_fix.sd_ix,2],
               fpr_fix[fpr_fix.sd_ix],
               tpr_quantile[fpr_fix.sd_ix,4])
      epsilon = 0.01
      segments(fpr_fix[fpr_fix.sd_ix]-epsilon,
               tpr_quantile[fpr_fix.sd_ix,2],
               fpr_fix[fpr_fix.sd_ix]+epsilon,
               tpr_quantile[fpr_fix.sd_ix,2])
      segments(fpr_fix[fpr_fix.sd_ix]-epsilon,
               tpr_quantile[fpr_fix.sd_ix,4],
               fpr_fix[fpr_fix.sd_ix]+epsilon,
               tpr_quantile[fpr_fix.sd_ix,4])
      
      return( list(cbind(tpr_quantile[,3],fpr_fix),cbind(fpr_fix.sd_ix,tpr_quantile[,c(2,4)])) )
      
      
    }
    
    
  }
  
  
  
  #by threshold
  if( type=="threshold" )
  {
    tpr_list = lapply( data,function(x){return(x[,1])} )
    tpr_mat = do.call(cbind,tpr_list)
    tpr_mean = rowMeans(tpr_mat)
    tpr_sd = apply(tpr_mat,1,sd)
    
    
    fpr_list = lapply(data,function(x){return(x[,2])} )
    fpr_mat = do.call(cbind,fpr_list)
    fpr_mean = rowMeans(fpr_mat)
    fpr_sd = apply(fpr_mat,1,sd)
    
    lg = length(tpr_sd)
    st = floor(lg/10)
    sd_ix = seq(st,lg,by=st)
    
    #draw average ROC curve
    plot(fpr_mean,tpr_mean,"l",
         xlim=c(0,1),ylim=c(0,1),
         xlab="False Postive Rate",
         ylab="True Positive Rate",
         main=title_str)
    lines(seq(0,1,0.1),seq(0,1,0.1),lty=2,col="gray")
    
    #vertical sd error bar
    segments(fpr_mean[sd_ix],tpr_mean[sd_ix]-tpr_sd[sd_ix],
             fpr_mean[sd_ix],tpr_mean[sd_ix]+tpr_sd[sd_ix])
    epsilon = 0.01
    segments(fpr_mean[sd_ix]-epsilon,tpr_mean[sd_ix]-tpr_sd[sd_ix],
             fpr_mean[sd_ix]+epsilon,tpr_mean[sd_ix]-tpr_sd[sd_ix])
    segments(fpr_mean[sd_ix]-epsilon,tpr_mean[sd_ix]+tpr_sd[sd_ix],
             fpr_mean[sd_ix]+epsilon,tpr_mean[sd_ix]+tpr_sd[sd_ix])
    #horiZontal sd error bar
    segments(fpr_mean[sd_ix]-fpr_sd[sd_ix],tpr_mean[sd_ix],
             fpr_mean[sd_ix]+fpr_sd[sd_ix],tpr_mean[sd_ix])
    segments(fpr_mean[sd_ix]-fpr_sd[sd_ix],tpr_mean[sd_ix]-epsilon,
             fpr_mean[sd_ix]-fpr_sd[sd_ix],tpr_mean[sd_ix]+epsilon)
    segments(fpr_mean[sd_ix]+fpr_sd[sd_ix],tpr_mean[sd_ix]-epsilon,
             fpr_mean[sd_ix]+fpr_sd[sd_ix],tpr_mean[sd_ix]+epsilon)
    
    return(list(cbind(tpr_mean,fpr_mean),cbind(tpr_sd,fpr_sd)))
  }
  
  
  
  
}


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

###load data###
##lap-top##
# pan_path = "C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/cnv/cisplatin/randomForest/pool/"
# info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2.5_fold_cv.mat.txt"
# curr_title = "Pan-cancer Analysis based on CNV" 
# pdf_file = "pan.cnv.rf.pdf"
#cluster#
args <- commandArgs(trailingOnly=TRUE)
pan_path = args[1]
info_file = args[2]
curr_title = args[3]
pdf_file = args[4]
##both##
cisplatin.info = read.table(info_file,header=T,quote="",sep="\t")
pan.files = list.files(path=pan_path,full.names=T,
                       pattern="pan.elanet.mid_res.test_[0-9]*.20150713.txt")
pan_score = lapply(pan.files,read_table)
filenames = lapply(pan.files,find_name)
names(pan_score) = filenames
###plot curves###
pdf(pdf_file)
##ROC curve##
fun_tp <- function(x){
  pats = rownames(x)
  x = as.vector(x[,1])
  names(x) = pats
  return(x) 
}
pan_score = lapply(pan_score,fun_tp)

cl = makeCluster(no_cores)
registerDoParallel(cl)
roc_res <- foreach(i=1:length(pan_score) ) %dopar%
{
  #for glmnet, second class is target class
  curr_pats = names(pan_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  #roc_res[[i]] = ensemble_roc(pan_score[[i]],resp,"sensitive")
  return( ensemble_roc(pan_score[[i]],resp,"insensitive") )
}
stopImplicitCluster()
stopCluster(cl)

ave_roc(roc_res,type="vertical",curr_title,"sd")
title(NULL,"Vertical average;std bar")
ave_roc(roc_res,type="vertical",curr_title,"quantile")
title(NULL,"Vertical average;quantile bar")
ave_roc(roc_res,type="threshold",curr_title,"Sd")
title(NULL,"Threshold average;std bar")

test_fit = auc_random(roc_res,test="ttest")
test_fit = auc_random(roc_res,test="ztest")
test_fit = auc_random(roc_res,test="wilcox")


##plot error rate##
cl = makeCluster(no_cores)
registerDoParallel(cl)
error_res <-foreach(i=1:length(pan_score)) %dopar%
{
  curr_pats = rownames(pan_score[[i]])
  resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
  #error_res[[i]] = identify_class(pan_score[[i]],resp,"sensitive")
  return( identify_class(pan_score[[i]],resp,"sensitive") )
}
stopImplicitCluster()
stopCluster(cl)
error_calc(error_res,cisplatin.info)

dev.off()
