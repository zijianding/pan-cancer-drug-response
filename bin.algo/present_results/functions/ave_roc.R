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