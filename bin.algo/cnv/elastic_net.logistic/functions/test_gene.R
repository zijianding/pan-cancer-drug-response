test_gene <- function(train_dat, test_dat, info, sig_gene = 1,
                      p_thresh=0.05,q_thresh=0.05, type="ttest",
                      p_step=0.001,q_step=0.05, p_up = 0.05,q_up = 0.2,
                      parallel=T,multi_cancer=T)
{
  #input data: row as genes and col as samples
  #output data:row as genes and col as samples
  #output significant threshold: q_thresh and p_thresh
  #type: default is t test; 
  #wilcox rank sum test, and logistic regression are also provided: "wilcox", "regression"
  #sig_gene, when no significant genes according to p/q value threshold
  
  
  if( type=="ttest" )
  {
    #perform on train data
    train_pats = colnames(train_dat)
    curr_info = info[match(train_pats,as.character(info$patient)),]
    sen_ix = which(as.character(curr_info$response)=="sensitive")
    insen_ix = which(as.character(curr_info$response)=="insensitive")
    
    #perform test
    p_values = vector( length=nrow(train_dat),mode="numeric" )
    if( parallel==F)
    {
      for( i in 1:nrow(train_dat) )
      {
        sen_dat = as.numeric(train_dat[i,sen_ix])
        insen_dat = as.numeric(train_dat[i,insen_ix])
        test_fit = t.test(sen_dat,insen_dat,na.action="na.omit")
        p_values[i] = test_fit$p.value
      }
    }
    if( parallel==T)
    {
      p_values <- foreach( i=1:nrow(train_dat),
                           .combine='c') %dopar%{
                             
                             sen_dat = as.numeric(train_dat[i,sen_ix])
                             insen_dat = as.numeric(train_dat[i,insen_ix])
                             test_fit = t.test(sen_dat,insen_dat,na.action="na.omit")
                             return(test_fit$p.value)
                           }
    }
    
    
    #multiple correction
    q_values = p.adjust(p_values,method="fdr",n=length(p_values))
    
    #determine q value threshold or p value threshold
    q_use = TRUE
    while( sum(q_values<q_thresh) <= sig_gene )
    {
      q_thresh = q_thresh + p_step
      if( q_thresh > q_up)
      {
        q_use = FALSE
        break
      }
    }
    q_use = FALSE #use p value at last
    p_use = TRUE
    if( q_use == FALSE )
    {
      while( sum(p_values<p_thresh) <= sig_gene  )
      {
        p_thresh = p_thresh + p_step
        if(p_thresh > p_up )
        {
          p_use = FALSE
          break
        }
      }
    }
    
    #identify significant genes
    if( q_use == TRUE )
    {
      sig_ix = q_values<=q_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"q_thresh",q_thresh))
    }
    if( (q_use==F) && (p_use==T) )
    {
      sig_ix = p_values <= p_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"p_thresh",p_thresh))
    }
    if( (q_use==F) && (p_use==F) )
    {
      sig_ix = sort(p_values,decreasing=F,index.return=T)$ix
      sig_ix = sig_ix[1:sig_gene]
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"sig_gene",sig_gene))
    }
    
  }
  if( type=="wilcox" )
  {
    #train data
    train_pats = colnames(train_dat)
    curr_info = info[match(train_pats,as.character(info$patient)),]
    sen_ix = which(as.character(curr_info$response)=="sensitive")
    insen_ix = which(as.character(curr_info$response)=="insensitive")
    
    #perform test
    p_values = vector( length=nrow(train_dat),mode="numeric" )
    if( parallel==F )
    {
      for( i in 1:nrow(train_dat) )
      {
        sen_dat = as.numeric(as.character(train_dat[i,sen_ix]))
        insen_dat = as.numeric(as.character(train_dat[i,insen_ix]))
        test_fit = wilcox.test(sen_dat,insen_dat,na.action="na.omit")
        p_values[i] = test_fit$p.value
      }
    }
    if( parallel==T)
    {
      p_values <- foreach( i=1:nrow(train.dat), .combine='c' ) %dopar%{
        sen_dat = as.numeric(as.character(train_dat[i,sen_ix]))
        insen_dat = as.numeric(as.character(train_dat[i,insen_ix]))
        test_fit = wilcox.test(sen_dat,insen_dat,na.action="na.omit")
        return(test_fit$p.value)
      }
    }
    
    #multiple correction
    q_values = p.adjust(p_values,method="fdr",n=length(p_values))
    
    #determine q value threshold or p value threshold
    q_use = TRUE
    while( sum(q_values<q_thresh) <= sig_gene )
    {
      q_thresh = q_thresh + q_step
      if( q_thresh > q_up)
      {
        q_use = FALSE
        break
      }
    }
    q_use = FALSE #use p value at last
    p_use = TRUE
    if( q_use == FALSE )
    {
      while( sum(p_values<p_thresh) <= sig_gene  )
      {
        p_thresh = p_thresh + p_step
        if(p_thresh > p_up )
        {
          p_use = FALSE
          break
        }
      }
    }
    
    #identify significant genes
    if( q_use == TRUE )
    {
      sig_ix = q_values<=q_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"q_thresh",q_thresh))
    }
    if( (q_use==F) && (p_use==T) )
    {
      sig_ix = p_values <= p_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"p_thresh",p_thresh))
    }
    if( (q_use==F) && (p_use==F) )
    {
      sig_ix = sort(p_values,decreasing=F,index.return=T)$ix
      sig_ix = sig_ix[1:sig_gene]
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"sig_gene",sig_gene))
    }
    
  }
  if( type =="regression")
  {
    #train data
    train_pats = colnames(train_dat)
    cancers = as.character(info$cancer[match(train_pats,as.character(info$patient))])
    cancers = as.factor(cancers)
    responses = as.character(info$response[match(train_pats,as.character(info$patient))])
    responses = as.factor(responses)
    
    #perform test
    p_values = vector( length=nrow(train_dat),mode="numeric" )
    if( parallel==F )
    {
      for(i in 1:nrow(train_dat) )
      {
        if(multi_cancer==T){
          fit = glm(responses~as.numeric(as.character(train_dat[i,]))+cancers,family=binomial)
        }else{
          fit = glm(responses~as.numeric(as.character(train_dat[i,])),family=binomial)
        }
        
        sum_fit = summary(fit)
        coeff = sum_fit$coefficients
        p_values[i] = coeff[2,ncol(coeff)]
      }
    }
    if( parallel==T )
    {
      p_values <- foreach(i=1:nrow(train_dat), .combine='c') %dopar%{
        if(multi_cancer==T){
          fit = glm(responses~as.numeric(as.character(train_dat[i,]))+cancers,family=binomial)
        }else{
          fit = glm(responses~as.numeric(as.character(train_dat[i,])),family=binomial)
        }
        sum_fit = summary(fit)
        coeff = sum_fit$coefficients
        return( coeff[2,ncol(coeff)] )
      }
      
    }
    
    #multiple correction
    q_values = p.adjust(p_values,method="fdr",n=length(p_values))
    
    #determine q value threshold or p value threshold
    q_use = TRUE
    while( sum(q_values<q_thresh) <= sig_gene )
    {
      q_thresh = q_thresh + q_step
      if( q_thresh > q_up )
      {
        q_use = FALSE
        break
      }
    }
    q_use = FALSE #use p value at last
    p_use = TRUE
    if( q_use == FALSE )
    {
      while( sum(p_values<p_thresh) <= sig_gene  )
      {
        p_thresh = p_thresh + p_step
        if(p_thresh > p_up)
        {
          p_use = FALSE
          break
        }
      }
    }
    
    #identify significant genes
    if( q_use == TRUE )
    {
      sig_ix = q_values<=q_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"q_thresh",q_thresh))
    }
    if( (q_use==F) && (p_use==T) )
    {
      sig_ix = p_values <= p_thresh
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"p_thresh",p_thresh))
    }
    if( (q_use==F) && (p_use==F) )
    {
      sig_ix = sort(p_values,decreasing=F,index.return=T)$ix
      sig_ix = sig_ix[1:sig_gene]
      train_dat = train_dat[sig_ix,]
      test_dat = test_dat[sig_ix,]
      return(list(train_dat,test_dat,"sig_gene",sig_gene))
    }
  }
  
  
  
  
}