####for mRNA-seq data###
###system parameters###
#options(digits = 15)

####load data####
#cluster
args <- commandArgs(trailingOnly=TRUE)
data_file = args[1]
info_file = args[2]
output_folder = args[3]#no "/" at the end
create_folder = args[4]

test_fold = as.numeric(as.character(args[5]))#the current test fold, ranging from 1 to 100
shuffle = as.numeric(as.character(args[6])) #"NULL"/1-100

input_type = args[7]   # molecular_only/clinical_only/clinical_molecular/***/***
output_type = args[8]  # "performance" for ROC, "marker" for extract marker, "shuffle" for permutation
calc_cancer = args[9]  # "sin_cancer"/"pan_cancer"
calc_gene = args[10]   # "all_gene"/"gene_set"

core.cancer = args[11] # "BLCA" etc/NULL, or train cancer for cross tumor
test.cancer = args[12] # test cancer for cross tumor
gene_set = args[13]    # NULL/gene_set/shuffle_gene_set(along with "shuffle")


source("source_all.R") #function file and parameter file

#desktop
# data_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2_focal.gdac_20141206.preprocess.txt"
# info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/cisplatin.gistic2_focal.5_fold_cv.mat.txt"
# output_folder = "C:/Users/zding/workspace/projects/drug_sensitivity/tmp/"
# create_folder = "test_cnv_map_gene"
# test_fold=1
# shuffle = NULL
# 
# input_type = "half_clinical_molecular" #NOTICE, input_type and output_type must be afront of source
# output_type = "performance"
# calc_cancer = "pan_cancer"
# calc_gene = "gene_set"
# 
# core.cancer = NULL
# test.cancer = NULL
# gene_set = "C:/Users/zding/workspace/projects/drug_sensitivity/data/gene_sets/cytoscape_cisplatin_gene.txt"
# 
# 
# source("source_all.R")



#both
cisplatin.dat = read.table(data_file,header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(info_file,sep="\t",header=T,quote="")
test_fold = test_fold + info_col




###preprocess data###
#core.info = cisplatin.info[as.character(cisplatin.info$cancer) %in% core.cancer,]
if( calc_cancer == "cross_cancer")
{
  train.cancer = core.cancer
  core.info.1 = cisplatin.info[as.character(cisplatin.info$cancer)==train.cancer,]
  core.info.2 = cisplatin.info[as.character(cisplatin.info$cancer)==test.cancer,]
  core.info = rbind(core.info.1,core.info.2)
  multi_cancer = F
}
if( calc_cancer == "sin_cancer")
{
  core.info = cisplatin.info[as.character(cisplatin.info$cancer)==core.cancer,]
}
if( calc_cancer == "pan_cancer")
{
  core.info = cisplatin.info
  #core.cancer = c("CESC","LUAD", "BLCA")
}
if( calc_gene == "gene_set" )
{
  pre_genes = read.table(gene_set,header=F,sep="\t",quote="")
  
  #cisplatin.dat = map_rna_gene(cisplatin.dat, pre_genes$V1)
  cisplatin.dat = map_cnv_gene(cisplatin.dat,pre_genes$V1)
}



##find data##
if( output_type!="marker"  ) # shuffle/performance
{
  if(calc_cancer=="cross_cancer"){
    #test data
    test.pats = as.character(core.info$patient[as.character(core.info$cancer)==test.cancer])
    test.info = core.info[as.character(core.info$cancer)==test.cancer,]
    test.dat = as.matrix(cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))])
    test.resp = as.character(test.info$response[match(test.pats,as.character(test.info$patient))])
    
    #train data
    train.pats = as.character(core.info$patient[as.character(core.info$cancer)==train.cancer])
    train.pats = sample(train.pats,size=length(train.pats),replace=FALSE)
    train.info = core.info[as.character(core.info$cancer)==train.cancer,]
    train.dat = as.matrix(cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))])
    
    if(output_type == "shuffle") #cross tumor permutation
    {
      #core.info.1; core.info; train.info
      random.resp = sample(as.character(core.info.1$response),replace=F)
      core.info.1$response = random.resp
      core.info = rbind(core.info.1, core.info.2)
      train.info = core.info.1
    }
  }else{
    #test data
    test.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="validation"])
    test.info = core.info[as.character(core.info[,test_fold])=="validation",]
    test.dat = as.matrix(cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))])
    test.resp = as.character(test.info$response[match(test.pats,as.character(test.info$patient))])
    
    #train data
    train.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="train"])
    train.pats = sample(train.pats,size=length(train.pats),replace=FALSE)
    train.info = core.info[as.character(core.info[,test_fold])=="train",]
    train.dat = as.matrix(cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))])
  }
  
  #delete all zero genes
  tmp_list = find_0(train.dat,test.dat)
  train.dat = tmp_list[[1]]
  test.dat = tmp_list[[2]]
  
  ##filter lowly expressed genes##
  if( filter_low_exp == TRUE )
  {
    data.tmp = filter_mRNA(cisplatin.dat, train.pats, test.pats, 
                           low.thresh="Q1", type="dispersion" )
    train.dat = data.tmp[[1]]; test.dat = data.tmp[[2]]
  }
  
  
  ##normalization in each cancer##
  if( exp_normalize == TRUE )
  {
    data.tmp = exp_norm(cbind(train.dat,test.dat), train.pats, test.pats, cisplatin.info)
    train.dat = data.tmp[[1]]
    test.dat = data.tmp[[2]]
  }
  
  
  ##select differential genes across cancer types##
  if( find_diff_genes == TRUE )
  {
    cl = makeCluster(no_cores)
    registerDoParallel(cl)
    list_tmp = test_gene(train.dat, test.dat, train.info,parallel=T,
                         type = test_type,sig_gene = sig_gene,
                         p_thresh=p_thresh,q_thresh=q_thresh,p_step=p_step,
                         q_step=q_step,p_up = p_up,q_up = q_up, multi_cancer=multi_cancer)
    stopImplicitCluster()
    stopCluster(cl)
    train.dat = list_tmp[[1]]
    test.dat = list_tmp[[2]]
    type = list_tmp[[3]]
    thresh = list_tmp[[4]]
    tmp_str = paste("With ",type," and threshold ",thresh,", ",
                    nrow(train.dat)," genes are remained",sep="")
    print(tmp_str)
    key_param[1,2] = thresh;
    key_param[2,2] = nrow(train.dat)
    
  }
  
  
  #impute NAs#
  tmp_list = impute_NA(train.dat,test.dat,nrow(train.dat))
  train.dat = tmp_list[[1]]
  test.dat = tmp_list[[2]]
  
  #add cancer type#
  if( add_clinic == TRUE )
  {
    tmp_list = cancer_dummy( train.dat, cisplatin.info)
    train.dat = tmp_list[[1]]
    dummy = tmp_list[[2]]
    test.dat = dummy_to_test(test.dat, cisplatin.info, dummy)
  }
}

if( output_type == "marker" )
{
  #test data
  test.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="validation"])
  test.info = core.info[as.character(core.info[,test_fold])=="validation",]
  test.dat = as.matrix(cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))])
  test.resp = as.character(test.info$response[match(test.pats,as.character(test.info$patient))])
  
  #train data
  train.pats = as.character(core.info$patient)
  train.pats = sample(train.pats,size=length(train.pats),replace=FALSE)
  train.info = core.info
  train.dat = as.matrix(cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))])
  
  #delete all zero genes
  tmp_list = find_0(train.dat,test.dat)
  train.dat = tmp_list[[1]]
  test.dat = tmp_list[[2]]
  
  ##filter lowly expressed genes##
  if( filter_low_exp == TRUE )
  {
    data.tmp = filter_mRNA(cisplatin.dat, train.pats, test.pats, 
                           low.thresh="Q1", type="dispersion" )
    train.dat = data.tmp[[1]]; test.dat = data.tmp[[2]]
  }
  
  
  ##normalization in each cancer##
  if( exp_normalize == TRUE )
  {
    data.tmp = exp_norm(cbind(train.dat,test.dat), train.pats, test.pats, cisplatin.info)
    train.dat = data.tmp[[1]]
    test.dat = data.tmp[[2]]
  }
  
  
  ##select differential genes across cancer types##
  if( find_diff_genes == TRUE )
  {
    cl = makeCluster(no_cores)
    registerDoParallel(cl)
    list_tmp = test_gene(train.dat, test.dat, cisplatin.info,parallel=T,
                         type = test_type,sig_gene = sig_gene,
                         p_thresh=p_thresh,q_thresh=q_thresh,p_step=p_step,
                         q_step=q_step,p_up = p_up,q_up = q_up, multi_cancer=multi_cancer)
    stopImplicitCluster()
    stopCluster(cl)
    train.dat = list_tmp[[1]]
    test.dat = list_tmp[[2]]
    type = list_tmp[[3]]
    thresh = nrow(train.dat)
    tmp_str = paste("With ",type," and threshold ",thresh,", ",
                    nrow(train.dat)," genes are remained",sep="")
    print(tmp_str)
  }
  
  
  #impute NAs#
  tmp_list = impute_NA(train.dat,test.dat,nrow(train.dat))
  train.dat = tmp_list[[1]]
  test.dat = tmp_list[[2]]
  
  #add cancer type#
  if( add_clinic == TRUE )
  {
    tmp_list = cancer_dummy( train.dat, cisplatin.info)
    train.dat = tmp_list[[1]]
    dummy = tmp_list[[2]]
    test.dat = dummy_to_test(test.dat, cisplatin.info, dummy)
  }
  
}

if( (input_type == "clinical_only") || (input_type=="clinical_pred") )
{
  #test data
  test.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="validation"])
  test.info = core.info[as.character(core.info[,test_fold])=="validation",]
  test.dat = as.matrix(cisplatin.dat[,match(test.pats,colnames(cisplatin.dat))])
  test.resp = as.character(test.info$response[match(test.pats,as.character(test.info$patient))])
  
  #train data
  train.pats = as.character(core.info$patient[as.character(core.info[,test_fold])=="train"])
  train.pats = sample(train.pats,size=length(train.pats),replace=FALSE)
  train.info = core.info[as.character(core.info[,test_fold])=="train",]
  train.dat = as.matrix(cisplatin.dat[,match(train.pats,colnames(cisplatin.dat))])
  
  #delete all zero genes
  tmp_list = find_0(train.dat,test.dat)
  train.dat = tmp_list[[1]]
  test.dat = tmp_list[[2]]
  
  #add cancer type and delete molecular data#
  gene_num = nrow(train.dat)
  tmp_list = cancer_dummy( train.dat, cisplatin.info)
  train.dat = tmp_list[[1]]
  dummy = tmp_list[[2]]
  test.dat = dummy_to_test(test.dat, cisplatin.info, dummy)
  
  train.dat = train.dat[-seq(1,gene_num,by=1),]
  test.dat = test.dat[-seq(1,gene_num,by=1),]
}


#train patients response
train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])

if( input_type != "clinical_pred")
{
  ###estimated robust features###
  #bootstrap samples#
  list.bs_mat = bootstrap_sample( colnames(train.dat), train.resp, BS=BS)
  bs_mat.pats = list.bs_mat[[1]]
  bs_mat.resp = list.bs_mat[[2]]
  
  ##estimates recurrent features##
  #model training#
  best_lambda = matrix(NA,nrow=length(alphas),ncol=BS)
  best_auc = matrix(NA,nrow=length(alphas),ncol=BS)
  for(i in 1:length(alphas))
  {
    cl = makeCluster(no_cores)
    registerDoParallel(cl)
    train_list <- foreach( bs=1:BS, .packages="glmnet",
                           .combine='comb', .multicombine=TRUE,
                           .init=list(list(), list()) ) %dopar%
    {
      #bootstrap sample data
      curr.train_dat = train.dat[ ,match(bs_mat.pats[,bs],colnames(train.dat)) ]
      curr.train_resp = bs_mat.resp[,bs]
      
      #delete all zero genes
      curr.test_dat = test.dat
      tmp_list = find_0(curr.train_dat,curr.test_dat)
      curr.train_dat = tmp_list[[1]]
      curr.test_dat = tmp_list[[2]]
      
      
      #record the best models in each bootstrap sample
      cv_fit = cv.glmnet( t(curr.train_dat), as.factor(curr.train_resp), 
                          family="binomial", type.measure="auc" )
      
      ix = match(cv_fit$lambda.1se,cv_fit$lambda)
      #beta = cv_fit$glmnet.fit$beta[,ix]
      auc = cv_fit$cvm[ix]
      return( list(cv_fit$lambda.1se,auc) )
    }
    stopImplicitCluster()
    stopCluster(cl)
    
    curr_lambda = unlist(train_list[[1]])
    curr_auc = unlist( train_list[[2]] )
    
    best_lambda[i,] = curr_lambda
    best_auc[i,] = curr_auc
    
  }
    
  #find best models and corresponding features
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  best_beta <- foreach(bs=1:BS,.combine='cbind',
                       .packages="glmnet") %dopar%
  {
    #best model
    ix = which.max(best_auc[,bs])
    alpha = alphas[ix]
    lambda = best_lambda[ix,bs]
    
    #bootstrap data
    curr.train_dat = train.dat[ ,match(bs_mat.pats[,bs],colnames(train.dat)) ]
    curr.train_resp = bs_mat.resp[,bs]
    
    #delete all zero genes
    curr.test_dat = test.dat
    tmp_list = find_0(curr.train_dat,curr.test_dat)
    curr.train_dat = tmp_list[[1]]
    curr.test_dat = tmp_list[[2]]
    
    #train model
    glm_fit = glmnet( t(curr.train_dat), as.factor(curr.train_resp), 
                      family="binomial",alpha=alpha,lambda=lambda)
    
    curr_beta = as.vector(glm_fit$beta)
    
    return(curr_beta)
  }
  stopImplicitCluster()
  stopCluster(cl)


  #recurrent features#
  if( select_gene == "by_freq")
  {
    rownames(best_beta) = rownames(train.dat)
    list_features = gene_selection(t(best_beta),freq=freq)
    while( length(list_features[[2]]) <= feature_min )
    {
      freq = freq - freq_step
      list_features = gene_selection(t(best_beta),freq=freq)
    }
    tmp_str = paste( "Frequency ",freq," to select",length(list_features[[2]]), "recurrent features" )
    print(tmp_str)
    key_param[3,2] = freq
    key_param[4,2] = length(list_features[[2]])
    
    #refine data by recurrent features
    train.dat = train.dat[list_features[[2]],]
    test.dat = test.dat[list_features[[2]],]
    feature_freq = list_features[[3]]
  }
  if( select_gene == "freqNweight")
  {
    rownames(best_beta) = rownames(train.dat)
    list_features = gene_selection_sd(best_beta,sd=sd)
    while( length(list_features[[2]]) < feature_min )
    {
      sd = sd - sd_step
      list_features = gene_selection_sd(best_beta,sd=sd)
    }
    tmp_str = paste("Sd",sd,"to select",length(list_features[[2]]),"recurrent features",sep=" ")
    print(tmp_str)
    key_param[3,2] = sd
    key_param[4,2] = length(list_features[[2]])
    
    #refine data by recurrent features
    train.dat = train.dat[list_features[[2]],]
    test.dat = test.dat[list_features[[2]],]
    feature_freq = list_features[[3]]
  }
  
  
}



##predict response by ensemble model
if( output_type == "marker" )
{
  #create output file folder#
  dir.create(file.path(output_folder, create_folder), showWarnings = FALSE)
  setwd(file.path(output_folder,create_folder))
  
  #selected features#
  tmp_str = "marker_molecular_only.txt"
  #write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
  if(select_gene == "freqNweight"){
    write.table(feature_freq,tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
  }else{
    write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
  }
  
  #key parameters#
  tmp_str = paste("pan.elanet.param.test_",test_fold-3,".20150701.txt",sep="")
  write.table(key_param,tmp_str,row.names=F,col.names=T,quote=F,sep="\t")
}
if( output_type != "marker" )
{
  if(combine_clinic==TRUE)
  {
    tmp_list = cancer_dummy( train.dat, cisplatin.info)
    train.dat = tmp_list[[1]]
    dummy = tmp_list[[2]]
    test.dat = dummy_to_test(test.dat, cisplatin.info, dummy)
  }
  
  
  ##refit models and test##
  #bootstrap samples#
  #train patients response
  train.resp = as.character(train.info$response[match(colnames(train.dat),as.character(train.info$patient))])
  list.bs_mat = bootstrap_sample( colnames(train.dat), train.resp, BS=BS)
  bs_mat.pats = list.bs_mat[[1]]
  bs_mat.resp = list.bs_mat[[2]]
  
  
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  test_list <- foreach(bs=1:BS,.packages="glmnet") %dopar%
  {
    #bootstrap sample
    curr.train_dat = train.dat[ ,match(bs_mat.pats[,bs],colnames(train.dat)) ]
    curr.train_resp = as.factor(bs_mat.resp[,bs])
    
    #delete all zero genes
    curr.test_dat = test.dat
    tmp_list = find_0(curr.train_dat,curr.test_dat)
    curr.train_dat = tmp_list[[1]]
    curr.test_dat = tmp_list[[2]]
    
    #fit the best model,should use the logistic regression
    glm_fit = glmnet(x=t(curr.train_dat),y=as.factor(curr.train_resp),
                     lambda=0, family="binomial")
    
    #predict on test data
    pred = predict( object=glm_fit,newx=t(test.dat), type="response" )
    
    return(as.vector(pred))
    
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  test_score = do.call(cbind,test_list)
  rownames(test_score) = test.pats
  
  
}
  

##output results##
if( output_type == "shuffle" )
{
  ###output results###
  #create output file folder#
  
  dir.create(file.path(output_folder, create_folder), showWarnings = FALSE)
  setwd(file.path(output_folder,create_folder))
  new_folder = paste(output_folder,create_folder,sep="/")
  dir.create(file.path(new_folder,as.character(shuffle)), showWarnings = F)
  setwd( file.path(new_folder,as.character(shuffle)) )
  
  #feature frequncy#
  tmp_str = paste("pan.elanet.feature_freq.test_",test_fold-3,".20150701.tiff",sep="")
  tiff(tmp_str)
  hist(feature_freq,main="Frequency of hittring for genes",xlab="hitting Freq",ylab="Freq",50)
  dev.off()
  
  ##final results##
  #selected features#
  tmp_str = paste("pan.elanet.feature.test_",test_fold-3,".20150701.txt",sep="")
  #write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
  if(select_gene == "freqNweight"){
    write.table(feature_freq,tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
  }else{
    write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
  }
  
  #test by each model#
  tmp_str = paste("pan.elanet.mid_res.test_",test_fold-3,".20150701.txt",sep="")
  write.table(t(test_score),tmp_str,col.names=T,row.names=F,sep="\t",quote=F)
  
  #key parameters#
  tmp_str = paste("pan.elanet.param.test_",test_fold-3,".20150701.txt",sep="")
  write.table(key_param,tmp_str,row.names=F,col.names=T,quote=F,sep="\t")
  
  ##plot TEST performance##
  tmp_str = paste("pan.elanat.test_",test_fold-3,".20150701.tiff",sep="")
  tiff(tmp_str)
  #roc
  roc = ensemble_roc(t(test_score),test.resp,"sensitive")
  plot(roc[,2],roc[,1],"b",cex=0.8,xlab="FPR",ylab="TPR",
       xlim=c(0,1),ylim=c(0,1))
  auc = trapz(roc[,2],roc[,1])
  tmp_str = paste("AUC = ",round(auc,digits=2),sep="")
  title("Test performance on Pan-Caner",tmp_str)
  dev.off()
  if(calc_cancer=="cross_cancer"){
    tmp_str = paste("pan.elanet.roc.test_",test_fold-3,".20150701.txt",sep="")
    colnames(roc) = c("tpr","fpr")
    write.table(roc,tmp_str,row.names=F,col.names=T,quote=F,sep="\t")
  }
}

if( output_type == "performance")
{
  ###output results###
  #create output file folder#
  dir.create(file.path(output_folder, create_folder), showWarnings = FALSE)
  setwd(file.path(output_folder,create_folder))
  
  
  
  ##final results##
  #selected features#
  if(input_type != "clinical_pred")
  {
    #feature frequncy#
    tmp_str = paste("pan.elanet.feature_freq.test_",test_fold-3,".20150701.tiff",sep="")
    tiff(tmp_str)
    hist(feature_freq,main="Frequency of hittring for genes",xlab="hitting Freq",ylab="Freq",50)
    dev.off()
    
    tmp_str = paste("pan.elanet.feature.test_",test_fold-3,".20150701.txt",sep="")
    #write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
    if(select_gene == "freqNweight"){
      write.table(feature_freq,tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
    }else{
      write.table(feature_freq[order(feature_freq,decreasing=T)],tmp_str,row.names=T,col.names=F,quote=F,sep="\t")
    }
  }
   
  
  #test by each model#
  tmp_str = paste("pan.elanet.mid_res.test_",test_fold-3,".20150701.txt",sep="")
  write.table(t(test_score),tmp_str,col.names=T,row.names=F,sep="\t",quote=F)
  
  #key parameters#
  tmp_str = paste("pan.elanet.param.test_",test_fold-3,".20150701.txt",sep="")
  write.table(key_param,tmp_str,row.names=F,col.names=T,quote=F,sep="\t")
  
  ##plot TEST performance##
  tmp_str = paste("pan.elanat.test_",test_fold-3,".20150701.tiff",sep="")
  tiff(tmp_str)
  #roc
  roc = ensemble_roc(t(test_score),test.resp,"sensitive")
  plot(roc[,2],roc[,1],"b",cex=0.8,xlab="FPR",ylab="TPR",
       xlim=c(0,1),ylim=c(0,1))
  auc = trapz(roc[,2],roc[,1])
  tmp_str = paste("AUC = ",round(auc,digits=2),sep="")
  title("Test performance on Pan-Caner",tmp_str)
  dev.off()
  if(calc_cancer=="cross_cancer"){
    tmp_str = paste("pan.elanet.roc.test_",test_fold-3,".20150701.txt",sep="")
    colnames(roc) = c("tpr","fpr")
    write.table(roc,tmp_str,row.names=F,col.names=T,quote=F,sep="\t")
  }
}


####THE END####






