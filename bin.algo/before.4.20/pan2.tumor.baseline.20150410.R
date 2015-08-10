#setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/methylation")
args <- commandArgs(trailingOnly=TRUE)
########################
######input data #######
###patient information###
#patient.info =  read.table("cancer.patient.drug.response.methylation.20150330.txt",header=T,sep="\t",quote="")
patient.info = read.table(args[1],header=T,sep="\t",quote="")
patient.info$response = as.character(patient.info$response)
patient.info = patient.info[as.character(patient.info$data)=="YES.methylation",]
tmp.arr = vector(mode="character",length=nrow(patient.info))
tmp.arr[which(as.character(patient.info$response)=="Clinical Progressive Disease")] = "insensitive"
tmp.arr[which(patient.info$response=="Stable Disease")] = "insensitive"
tmp.arr[which(patient.info$response=="Complete Response")] = "sensitive"
tmp.arr[which(patient.info$response=="Partial Response")] = "sensitive"
patient.info$response = tmp.arr
cisplatin.info = patient.info[as.character(patient.info$drug)=="Cisplatin",]
################
###input methylation data of the drug###
#patient.data = read.table("Cisplatin.methylation450k.gdac_20141206.preprocess.txt",sep="\t",header=T,quote="",row.names=1)
patient.data = read.table(args[2],sep="\t",header=T,quote="",row.names=1)
#delte the samples NOT primary tumor
primary.tumor.arr = c()
for( i in 2:ncol(patient.data) )
{
  tmp = strsplit(as.character(colnames(patient.data)[i]),"\\.")
  #tmp.cancer = tmp[[1]][1]
  #tmp.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
  #if( tmp[[1]][5] == "01A" )
  if( !is.na( pmatch("01",tmp[[1]][5]) ) )
  {
    primary.tumor.arr = c(primary.tumor.arr,i)
    #patients = c(patients,curr.pat)
  }
  else
  {
    print(tmp[[1]][5])
  }
}
patient.data = patient.data[,primary.tumor.arr] #then no column of gene
#patient.data = patient.data[,colnames(patient.data)]
#match the data and cisplatin.info
patients = c()
for( i in 1:ncol(patient.data) )
{
  tmp = strsplit(as.character(colnames(patient.data)[i]),"\\.")
  #cancer = tmp[[1]][1]
  tmp.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
  patients = c(patients,tmp.pat)
  #   #if( tmp[[1]][5] == "01A" )
  #   if( !is.na( pmatch("01",tmp[[1]][5]) ) )
  #   {
  #primary.tumor.arr = c(primary.tumor.arr,i)
  
  #   }
  #   else
  #   {
  #     print(tmp[[1]][5])
  #   }
}
cisplatin.info = cisplatin.info[match(patients,as.character(cisplatin.info$patient)),]


#each sample sensitive or not,index
# sen.index = c()
# insen.index = c()
# patients = colnames(patient.data)
# for(i in 1:length(patients))
# {
#   tmp = strsplit(patients[i],"\\.")
#   cancer = tmp[[1]][1]
#   curr.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
#   pat.ix = match(curr.pat,as.character(cisplatin.info$patient))
#   
#   if( cisplatin.info$response[pat.ix] == "sensitive")
#   {
#     sen.index = c(sen.index,i )
#   }
#   if( cisplatin.info$response[pat.ix] == "insensitive")
#   {
#     insen.index = c(insen.index,i)
#   }
# }
#core cancer
cancer.table = table(as.character(cisplatin.info$cancer))
cancer.table = cancer.table[order(cancer.table,decreasing=T)]
sen.arr = c()
insen.arr = c()
for( i in 1:length(cancer.table) )
{
  curr.cancer = names(cancer.table)[i]
  curr.df = cisplatin.info[as.character(cisplatin.info$cancer)==curr.cancer,]
  sen.num = sum(as.character(curr.df$response)=="sensitive")
  insen.num = sum(as.character(curr.df$response)=="insensitive")
  sen.arr = c(sen.arr,sen.num)
  insen.arr = c(insen.arr,insen.num)
}
cancer.df = data.frame(cancer.table,sen.arr,insen.arr)
colnames(cancer.df) = c("patient","sensitive","insensitive")
cancer.df = cancer.df[cancer.df[,1]!=1,]

#library(ggplot2)
#ggplot(cancer.df,aes(x=colnames(cancer.df),y=patient)) + geom_bar(aes(fill=sensitive),stat="identity")
#core cancer: LUAD, BLCA, MESO, STAD

#########################
####model constrution####
core.cancer = c("LUAD","BLCA","MESO","STAD")
library(glmnet)

#pick up the data of core cancer
core.ix = c()
for(i in 1:ncol(patient.data))
{
  tmp = strsplit(colnames(patient.data)[i],"\\.")
  cancer = tmp[[1]][1]
  if( !is.na(match(cancer,core.cancer)) )
  {
    core.ix = c(core.ix,i)
  }
}
core.ix = sample(core.ix,size=length(core.ix),replace=F)
core.dat = patient.data[,core.ix]
core.response = vector(length=ncol(core.dat),mode="character")
for(i in 1:ncol(core.dat))
{
  tmp = strsplit(colnames(core.dat)[i],"\\.")
  pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
  tmp.ix = match(pat,as.character(patient.info$patient))
  core.response[i] = as.character(patient.info$response[tmp.ix])
}

#10-fold CV
fd.size = round(ncol(core.dat)/10)
fd.ix = rep(1:10,each = fd.size)
tmp = ncol(core.dat)%%10
fd.ix = c(fd.ix,rep(10,tmp))
for( f in 1:10 )
{
  ##find train data
  train.ix = which(fd.ix!=f)
  test.ix = which(fd.ix==f)
  aov.pvalue = vector(length=nrow(core.dat),mode="numeric")
  probe.cancer = c()
  probe.resp = c()
  probe.ix =c()
  for(j in 1:ncol(core.dat))
  {
    tmp = strsplit(colnames(core.dat)[j],"\\.")
    #NOTICE Only train data
    if( !is.na(match(j,train.ix)))
    {
      probe.cancer = c(probe.cancer, tmp[[1]][1])
      probe.resp = c(probe.resp,core.response[j])
      #probe.val = c(probe.val,as.numeric(as.character(core.dat[i,j])))
      probe.ix = c(probe.ix,j)
    }
  }
  for(i in 1:nrow(core.dat))
  {
    probe.val = unlist(core.dat[i,probe.ix])
    #probe.df = data.frame(probe.val,probe.cancer,probe.resp)
    aov.res = summary(aov(probe.val~probe.resp*probe.cancer))
    aov.pvalue[i] =aov.res[[1]]$Pr[1]
  }
  aov.fdr = p.adjust(aov.pvalue,method="fdr")
  
  #if no significant probes
  if( sum(aov.fdr<=0.05)<=1 )
  {
    print(paste("NO sig probes in fold ",f,"of core cancer set",sep=" "))
    #the 1st 1000:10000,pace=10000 probes
    fdr.arr.ix = sort(aov.fdr,index.return=T)$ix
    probe.num = seq(1000,10000,by=1000)
    for( p in 1:length(probe.num) )
    {
      tmp.char = paste("The first",probe.num[p],"probes","of test fold",f,sep=" ")
      print(tmp.char)
      ########################################
      pre.probe.ix = fdr.arr.ix[1:probe.num[p]]
      #now.dat = core.dat[pre.probe.ix,]
      ###############################
      #training data and testing data
      train.dat = t(core.dat[pre.probe.ix,train.ix])
      for(i in 1:ncol(train.dat))
      {
        tmp = which(is.na(train.dat[,i]))
        if( length(tmp)>0 )
        {
          curr.mean = mean(train.dat[,i],na.rm=T)
          train.dat[tmp,i] = curr.mean
        }
      }
      train.response = core.response[train.ix]
      
      ##############################
      #construct logistic regression
      cv.res = cv.glmnet(x=train.dat,y=as.factor(train.response),alpha=1,nfolds=5,type.measure="class",family="binomial")
      plot(cv.res)
      
      glm.res = glmnet(x=train.dat,y=as.factor(train.response),lambda=cv.res$lambda.min,family="binomial",alpha=1)
      
      #####
      #prediction on test set
      test.dat = t(core.dat[pre.probe.ix,test.ix])
      for(i in 1:ncol(test.dat))
      {
        tmp = which(is.na(test.dat[,i]))
        if( length(tmp)>0 )
        {
          curr.mean = mean(test.dat[,i],na.rm=T)
          test.dat[tmp,i] = curr.mean
        }
      }
      test.response = core.response[test.ix]
      
      
      prediction.res = predict(object = glm.res,newx = test.dat,s=cv.res$lambda.min,type="response")
      prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
      prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(train.response))[2]
      prediction.res.arr[which(prediction.res<0.5)] = levels(as.factor(train.response))[1]
      
      
      #mis.rate = sum(prediction.res.arr!=test.response)
      #mis.rate = mis.rate/length(prediction.res.arr)
      ###the proprotion
      insen.ratio = sum(test.response=="insensitive")/length(test.response)
      tmp.char = paste("The proportion of insensitive samples of test fold",f,"is",insen.ratio,sep=" ")
      print(tmp.char)
      
      ###the predicted results
      tp=0
      fp=0
      fn=0
      tn=0
      if( sum(prediction.res.arr=="insensitive")>0 )
      {
        fp = sum(test.response[which(prediction.res.arr=="insensitive")]=="sensitive")
        tp = sum(test.response[which(prediction.res.arr=="insensitive")]=="insensitive")
      }
      if( sum(prediction.res.arr=="sensitive")>0 )
      {
        fn = sum(test.response[which(prediction.res.arr=="sensitive")]=="insensitive")
        tn = sum(test.response[which(prediction.res.arr=="sensitive")]=="sensitive")
      }
      #FDR = FP/(FP+TP), as small as possible
      if((tp+fp)>0)
      {
        fdr = fp/(fp+tp)
        tmp.char = paste("The FDR(small)",fdr,"of fold",f,"in core cancer set",sep=" ")
        print(tmp.char)
      }
      if((tp+fp)==0)
      {
        tmp.char = paste("The FDR(small) is NA of fold",f,"in core cancer set",sep=" ")
        print(tmp.char)
      }
      #Recall = TP/(TP+FN), as large as possible
      if((tp+fn)>0)
      {
        tmp.char = paste("The recall(large) is",tp/(tp+fn),"of fold",f,"in core cancer set",sep=" ")
        print(tmp.char)
      }
      if((tp+fn)==0)
      {
        tmp.char = paste("The recall(large) is NA of fold",f,"in core cancer set",sep=" ")
        print(tmp.char)
      }
      #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
      if((fp+tn)>0)
      {
        tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"of fold",f,"in core cancer set",sep=" ")
        print(tmp.char)
      }
      if((fp+tn)==0)
      {
        tmp.char = paste("The FPR(small) is NA of fold",f,"in core cancer set",sep=" ")
        print(tmp.char)
      }
      ###the misclassification rate
      mis.rate = (fp+fn)/(fp+fn+tn+tp)
      tmp.char = paste("The misclassification rate on test fold",f,"is",mis.rate,sep=" ")
      print(tmp.char)##the error, FDR, and recall for single tumor
      
      
      ###on core cancer, each as test
      test.core.cancer = vector(length=nrow(test.dat),mode="character")
      test.core.resp = vector(length=nrow(test.dat),mode="character")
      for(i in 1:nrow(test.dat))
      {
        tmp = strsplit(rownames(test.dat)[i],"\\.")
        test.core.cancer[i] = tmp[[1]][1]
        #tmp.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
        #tmp.ix = match(tmp.pat,as.character(patient.info$patient))
        #test.core.re
      }
      for(i in 1:length(core.cancer))
      {
        tmp.char = paste("Test cancer:",core.cancer[i],"by core cancer set of test fold",f,sep=" ")
        print(tmp.char)
        tmp.ix = which(test.core.cancer==core.cancer[i])
        pred = prediction.res.arr[tmp.ix]
        real = test.response[tmp.ix]
        
        insen.ratio = sum(real=="insensitive")/length(real)
        tmp.char = paste("The proportion of positive samples of test core",core.cancer[i] ,"of core cancer set is", insen.ratio,sep=" ")
        print(tmp.char)
        
        ###the predicted results
        tp=0
        fp=0
        fn=0
        tn=0
        if( sum(pred=="insensitive")>0 )
        {
          fp = sum(real[which(pred=="insensitive")]=="sensitive")
          tp = sum(real[which(pred=="insensitive")]=="insensitive")
        }
        if( sum(pred=="sensitive")>0 )
        {
          fn = sum(real[which(pred=="sensitive")]=="insensitive")
          tn = sum(real[which(pred=="sensitive")]=="sensitive")
        }
        #FDR = FP/(FP+TP), as small as possible
        if((tp+fp)>0)
        {
          fdr = fp/(fp+tp)
          tmp.char = paste("The FDR(small)",fdr,"on test core cancer",core.cancer[i],"set by test of fold",f,sep=" ")
          print(tmp.char)
        }
        if((tp+fp)==0)
        {
          tmp.char = paste("The FDR(small) is NA on test core cancer",core.cancer[i], "by test of fold",f,sep=" ")
          print(tmp.char)
        }
        #Recall = TP/(TP+FN), as large as possible
        if((tp+fn)>0)
        {
          tmp.char = paste("The recall(large) is",tp/(tp+fn),"on test core cancer",core.cancer[i],"by test of fold",f,sep=" ")
          print(tmp.char)
        }
        if((tp+fn)==0)
        {
          tmp.char = paste("The recall(large) is NA","on test core cancer",core.cancer[i], "by test of fold",f,sep=" ")
          print(tmp.char)
        }
        #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
        if((fp+tn)>0)
        {
          tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"on test core cancer",core.cancer[i],"by test of fold",f,sep=" ")
          print(tmp.char)
        }
        if((fp+tn)==0)
        {
          tmp.char = paste("The FPR(small) is NA on test core cancer",core.cancer[i], "by test of fold", f,sep=" ")
          print(tmp.char)
        }
        ###the misclassification rate
        mis.rate = (fp+fn)/(fp+fn+tn+tp)
        tmp.char = paste("The misclassification rate is",mis.rate,"on test core cancer",core.cancer[i],"by test of fold",f,sep=" ")
        print(tmp.char)
      }
      
      
      
      ###############
      ##on other tumor
      for(c in 1:nrow(cancer.df))
      {
        test.cancer = rownames(cancer.df)[c]
        #other tumor
        if( is.na(match(test.cancer,core.cancer)) )
        {
          tmp.char = paste("Test cancer:",test.cancer,"by core cancer set of test fold",f,sep=" ")
          print(tmp.char)
          #find samples for current test cancer
          #patients = colnames(patient.data)
          test.sen = c()
          test.insen = c()
          pat.ix = c()
          for(i in 1:ncol(patient.data) )
          {
            tmp = strsplit(colnames(patient.data)[i],"\\.")
            cancer = tmp[[1]][1]
            test.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
            tmp.ix = match(test.pat,as.character(cisplatin.info$patient))
            if( cancer == test.cancer )
            {
              pat.ix = c(pat.ix,i)
              response = as.character(cisplatin.info$response[tmp.ix])
              if( response == "sensitive" )
              {
                test.sen = c(test.sen,length(pat.ix))
              }
              if( response == "insensitive" )
              {
                test.insen = c( test.insen, length(pat.ix) )
              }
            }
          }
          test.response = vector(length=length(pat.ix),mode="character")
          test.dat = core.dat[pre.probe.ix,pat.ix]
          for(i in 1:ncol(test.dat))
          {
            tmp = strsplit(colnames(test.dat)[i],"\\.")
            #cancer = tmp[[1]][1]
            test.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
            tmp.ix = match(test.pat,patient.info$patient)
            test.response[i] = as.character(patient.info$response[tmp.ix])
          }
          test.response[test.sen] = "sensitive"
          test.response[test.insen] = "insensitive"
          
          #prediction part
          #find the probes
          test.dat = t(test.dat)
          for(i in 1:ncol(test.dat))
          {
            tmp = which(is.na(test.dat[,i]))
            if( length(tmp)>0 )
            {
              curr.mean = mean(test.dat[,i],na.rm=T)
              test.dat[tmp,i] = curr.mean
            }
          }
          
          prediction.res = predict(object = glm.res,newx = test.dat,s=cv.res$lambda.min,type="response")
          prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
          prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(train.response))[2]
          prediction.res.arr[which(prediction.res<0.5)] = levels(as.factor(train.response))[1]
          
          
          #mis.rate = sum(prediction.res.arr!=test.response)
          #mis.rate = mis.rate/length(prediction.res.arr)
          ###the proprotion
          insen.ratio = sum(test.response=="insensitive")/length(test.response)
          tmp.char = paste("The proportion of positive samples of test",test.cancer,"by fold",f,"is", insen.ratio,sep=" ")
          print(tmp.char)
          
          tp=0
          fp=0
          fn=0
          tn=0
          if( sum(prediction.res.arr=="insensitive")>0 )
          {
            fp = sum(test.response[which(prediction.res.arr=="insensitive")]=="sensitive")
            tp = sum(test.response[which(prediction.res.arr=="insensitive")]=="insensitive")
          }
          if( sum(prediction.res.arr=="sensitive")>0 )
          {
            fn = sum(test.response[which(prediction.res.arr=="sensitive")]=="insensitive")
            tn = sum(test.response[which(prediction.res.arr=="sensitive")]=="sensitive")
          }
          #FDR = FP/(FP+TP), as small as possible
          if((tp+fp)>0)
          {
            fdr = fp/(fp+tp)
            tmp.char = paste("The FDR(small)",fdr,"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
            print(tmp.char)
          }
          if((tp+fp)==0)
          {
            tmp.char = paste("The FDR(small) is NA","on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
            print(tmp.char)
          }
          #Recall = TP/(TP+FN), as large as possible
          if((tp+fn)>0)
          {
            tmp.char = paste("The recall(large) is",tp/(tp+fn),"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
            print(tmp.char)
          }
          if((tp+fn)==0)
          {
            tmp.char = paste("The recall(large) is NA","on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
            print(tmp.char)
          }
          #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
          if((fp+tn)>0)
          {
            tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
            print(tmp.char)
          }
          if((fp+tn)==0)
          {
            tmp.char = paste("The FPR(small) is NA","on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
            print(tmp.char)
          }
          ###the misclassification rate
          mis.rate = (fp+fn)/(fp+fn+tn+tp)
          tmp.char = paste("The misclassification rate is",mis.rate,"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
          print(tmp.char)
          
          #construct SVM
          # library(e1071)
          # svm.linear.res = svm(train.dat,train.response,type="C-classification",kernel="linear",cost=10)
          # svm.rbf.res = svm(train.dat,train.response,type="C-classification",kernel="radius",gamma=1,cost=10)
        }
      }
    } 
  }
  #if some signifcant probes
  else
  {
    sig.num = sum(aov.fdr<=0.05)
    tmp.char = paste(sig.num,"sig probes by core cancer set of fold",f,sep=" ")
    print(tmp.char)
    
    pre.probe.ix = which(aov.fdr<=0.05)
    train.dat = t(core.dat[pre.probe.ix,train.ix])
    train.response = core.response[train.ix]
    for(i in 1:ncol(train.dat))
    {
      tmp = which(is.na(train.dat[,i]))
      if( length(tmp)>0 )
      {
        curr.mean = mean(train.dat[,i],na.rm=T)
        train.dat[tmp,i] = curr.mean
      }
    }
    
    #train model
    cv.res = cv.glmnet(x=train.dat,y=as.factor(train.response),alpha=1,nfolds=5,type.measure="class",family="binomial")
    plot(cv.res)
    
    glm.res = glmnet(x=train.dat,y=as.factor(train.response),lambda=cv.res$lambda.min,family="binomial",alpha=1)
    
    #test model
    test.dat = t(core.dat[pre.probe.ix,test.ix])
    test.response = core.response[test.ix]
    for(i in 1:ncol(test.dat))
    {
      tmp = which(is.na(test.dat[,i]))
      if( length(tmp)>0 )
      {
        curr.mean = mean(test.dat[,i],na.rm=T)
        test.dat[tmp,i] = curr.mean
      }
    }
    
    prediction.res = predict(object = glm.res,newx = test.dat,s=cv.res$lambda.min,type="response")
    prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
    prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(train.response))[2]
    prediction.res.arr[which(prediction.res<0.5)] = levels(as.factor(train.response))[1]
    
    
    #mis.rate = sum(prediction.res.arr!=test.response)
    #mis.rate = mis.rate/length(prediction.res.arr)
    ###the proprotion
    insen.ratio = sum(test.response=="insensitive")/length(test.response)
    tmp.char = paste("The proportion of positive samples of test set of core cancer set is", insen.ratio,sep=" ")
    print(tmp.char)
    
    ###the predicted results
    tp=0
    fp=0
    fn=0
    tn=0
    if( sum(prediction.res.arr=="insensitive")>0 )
    {
      fp = sum(test.response[which(prediction.res.arr=="insensitive")]=="sensitive")
      tp = sum(test.response[which(prediction.res.arr=="insensitive")]=="insensitive")
    }
    if( sum(prediction.res.arr=="sensitive")>0 )
    {
      fn = sum(test.response[which(prediction.res.arr=="sensitive")]=="insensitive")
      tn = sum(test.response[which(prediction.res.arr=="sensitive")]=="sensitive")
    }
    #FDR = FP/(FP+TP), as small as possible
    if((tp+fp)>0)
    {
      fdr = fp/(fp+tp)
      tmp.char = paste("The FDR(small)",fdr,"on test core cancer set by test of fold",f,sep=" ")
      print(tmp.char)
    }
    if((tp+fp)==0)
    {
      tmp.char = paste("The FDR(small) is NA on test core cancer set by test of fold",f,sep=" ")
      print(tmp.char)
    }
    #Recall = TP/(TP+FN), as large as possible
    if((tp+fn)>0)
    {
      tmp.char = paste("The recall(large) is",tp/(tp+fn),"on test core cancer set by test of fold",f,sep=" ")
      print(tmp.char)
    }
    if((tp+fn)==0)
    {
      tmp.char = paste("The recall(large) is NA","on test core cancer set by test of fold",f,sep=" ")
      print(tmp.char)
    }
    #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
    if((fp+tn)>0)
    {
      tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"on test core cancer set by test of fold",f,sep=" ")
      print(tmp.char)
    }
    if((fp+tn)==0)
    {
      tmp.char = paste("The FPR(small) is NA","on test core cancer set by test of fold",f,sep=" ")
      print(tmp.char)
    }
    ###the misclassification rate
    mis.rate = (fp+fn)/(fp+fn+tn+tp)
    tmp.char = paste("The misclassification rate is",mis.rate,"on test core cancer set by test of fold",f,sep=" ")
    print(tmp.char)
    
    ##the error, FDR, and recall for single tumor
    test.core.cancer = vector(length=nrow(test.dat),mode="character")
    test.core.resp = vector(length=nrow(test.dat),mode="character")
    for(i in 1:nrow(test.dat))
    {
      tmp = strsplit(rownames(test.dat)[i],"\\.")
      test.core.cancer[i] = tmp[[1]][1]
      #tmp.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
      #tmp.ix = match(tmp.pat,as.character(patient.info$patient))
      #test.core.re
    }
    for(i in 1:length(core.cancer))
    {
      tmp.char = paste("Test cancer:",core.cancer[i],"by core cancer set of test fold",f,sep=" ")
      print(tmp.char)
      tmp.ix = which(test.core.cancer==core.cancer[i])
      pred = prediction.res.arr[tmp.ix]
      real = test.response[tmp.ix]
      
      insen.ratio = sum(real=="insensitive")/length(real)
      tmp.char = paste("The proportion of positive samples of test core",core.cancer[i] ,"of core cancer set is", insen.ratio,sep=" ")
      print(tmp.char)
      
      ###the predicted results
      tp=0
      fp=0
      fn=0
      tn=0
      if( sum(pred=="insensitive")>0 )
      {
        fp = sum(real[which(pred=="insensitive")]=="sensitive")
        tp = sum(real[which(pred=="insensitive")]=="insensitive")
      }
      if( sum(pred=="sensitive")>0 )
      {
        fn = sum(real[which(pred=="sensitive")]=="insensitive")
        tn = sum(real[which(pred=="sensitive")]=="sensitive")
      }
      #FDR = FP/(FP+TP), as small as possible
      if((tp+fp)>0)
      {
        fdr = fp/(fp+tp)
        tmp.char = paste("The FDR(small)",fdr,"on test core cancer",core.cancer[i],"set by test of fold",f,sep=" ")
        print(tmp.char)
      }
      if((tp+fp)==0)
      {
        tmp.char = paste("The FDR(small) is NA on test core cancer",core.cancer[i], "by test of fold",f,sep=" ")
        print(tmp.char)
      }
      #Recall = TP/(TP+FN), as large as possible
      if((tp+fn)>0)
      {
        tmp.char = paste("The recall(large) is",tp/(tp+fn),"on test core cancer",core.cancer[i],"by test of fold",f,sep=" ")
        print(tmp.char)
      }
      if((tp+fn)==0)
      {
        tmp.char = paste("The recall(large) is NA","on test core cancer",core.cancer[i], "by test of fold",f,sep=" ")
        print(tmp.char)
      }
      #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
      if((fp+tn)>0)
      {
        tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"on test core cancer",core.cancer[i],"by test of fold",f,sep=" ")
        print(tmp.char)
      }
      if((fp+tn)==0)
      {
        tmp.char = paste("The FPR(small) is NA on test core cancer",core.cancer[i], "by test of fold", f,sep=" ")
        print(tmp.char)
      }
      ###the misclassification rate
      mis.rate = (fp+fn)/(fp+fn+tn+tp)
      tmp.char = paste("The misclassification rate is",mis.rate,"on test core cancer",core.cancer[i],"by test of fold",f,sep=" ")
      print(tmp.char)
    }
    
    
    #construct SVM
    # library(e1071)
    # svm.linear.res = svm(train.dat,train.response,type="C-classification",kernel="linear",cost=10)
    # svm.rbf.res = svm(train.dat,train.response,type="C-classification",kernel="radius",gamma=1,cost=10)
    
    ############################
    #######on other tumor#######
    for(c in 1:nrow(cancer.df))
    {
      test.cancer = rownames(cancer.df)[c]
      #other tumor
      if( is.na(match(test.cancer,core.cancer)) )
      {
        tmp.char = paste("Test cancer:",test.cancer,"by core cancer set","of test fold",f,sep=" ")
        print(tmp.char)
        ####################################
        #find samples for current test cancer
        #patients = colnames(patient.data)
        test.sen = c()
        test.insen = c()
        pat.ix = c()
        for(i in 1:ncol(patient.data) )
        {
          tmp = strsplit(colnames(patient.data)[i],"\\.")
          cancer = tmp[[1]][1]
          test.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
          tmp.ix = match(test.pat,as.character(cisplatin.info$patient))
          if( cancer == test.cancer )
          {
            pat.ix = c(pat.ix,i)
            response = as.character(cisplatin.info$response[tmp.ix])
            if( response == "sensitive" )
            {
              test.sen = c(test.sen,length(pat.ix))
            }
            if( response == "insensitive" )
            {
              test.insen = c( test.insen, length(pat.ix) )
            }
          }
        }
        test.dat = patient.data[,pat.ix]
        test.response = vector(length=length(pat.ix),mode="character")
        if(length(pat.ix)==1)
        {
          tmp = strsplit(colnames(patient.data)[pat.ix],"\\.")
          #cancer = tmp[[1]][1]
          test.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
          tmp.ix = match(test.pat,patient.info$patient)
          test.response[i] = as.character(patient.info$response[tmp.ix])
          
          if( length(test.sen)>0)
          {
            test.response[test.sen] = "sensitive"
          }
          if( length(test.insen)>0 )
          {
            test.response[test.insen] = "insensitive"
          }
          
          test.dat = test.dat[pre.probe.ix]
          tmp = which(is.na(test.dat))
          if( length(tmp)>0 )
          {
            curr.mean = mean(test.dat,na.rm=T)
            test.dat[tmp] = curr.mean
          } 
          
        }
        if( length(pat.ix)>1)
        {
          for(i in 1:ncol(test.dat))
          {
            tmp = strsplit(colnames(test.dat)[i],"\\.")
            #cancer = tmp[[1]][1]
            test.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
            tmp.ix = match(test.pat,patient.info$patient)
            test.response[i] = as.character(patient.info$response[tmp.ix])
          }
          test.response[test.sen] = "sensitive"
          test.response[test.insen] = "insensitive"
          
          test.dat = t(test.dat[pre.probe.ix,])
          for(i in 1:ncol(test.dat))
          {
            tmp = which(is.na(test.dat[,i]))
            if( length(tmp)>0 )
            {
              curr.mean = mean(test.dat[,i],na.rm=T)
              test.dat[tmp,i] = curr.mean
            }
          }
        }
        
        
        #prediction part
        #find the probes
        
        
        prediction.res = predict(object = glm.res,newx = test.dat,s=cv.res$lambda.min,type="response")
        prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
        prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(train.response))[2]
        prediction.res.arr[which(prediction.res<0.5)] = levels(as.factor(train.response))[1]
        
        
        #mis.rate = sum(prediction.res.arr!=test.response)
        #mis.rate = mis.rate/length(prediction.res.arr)
        ###the proprotion
        insen.ratio = sum(test.response=="insensitive")/length(test.response)
        tmp.char = paste("The proportion of positive samples of test",test.cancer,"by test fold of",f,"on core cnacer set is", insen.ratio,sep=" ")
        print(tmp.char)
        
        ###the predicted results
        tp=0
        fp=0
        fn=0
        tn=0
        if( sum(prediction.res.arr=="insensitive")>0 )
        {
          fp = sum(test.response[which(prediction.res.arr=="insensitive")]=="sensitive")
          tp = sum(test.response[which(prediction.res.arr=="insensitive")]=="insensitive")
        }
        if( sum(prediction.res.arr=="sensitive")>0 )
        {
          fn = sum(test.response[which(prediction.res.arr=="sensitive")]=="insensitive")
          tn = sum(test.response[which(prediction.res.arr=="sensitive")]=="sensitive")
        }
        #FDR = FP/(FP+TP), as small as possible
        if((tp+fp)>0)
        {
          fdr = fp/(fp+tp)
          tmp.char = paste("The FDR(small)",fdr,"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
          print(tmp.char)
        }
        if((tp+fp)==0)
        {
          tmp.char = paste("The FDR(small) is NA","on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
          print(tmp.char)
        }
        #Recall = TP/(TP+FN), as large as possible
        if((tp+fn)>0)
        {
          tmp.char = paste("The recall(large) is",tp/(tp+fn),"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
          print(tmp.char)
        }
        if((tp+fn)==0)
        {
          tmp.char = paste("The recall(large) is NA","on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
          print(tmp.char)
        }
        #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
        if((fp+tn)>0)
        {
          tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
          print(tmp.char)
        }
        if((fp+tn)==0)
        {
          tmp.char = paste("The FPR(small) is NA","on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
          print(tmp.char)
        }
        ###the misclassification rate
        mis.rate = (fp+fn)/(fp+fn+tn+tp)
        tmp.char = paste("The misclassification rate is",mis.rate,"on",test.cancer,"by test fold of",f,"core cancerset",sep=" ")
        print(tmp.char)
        
        #construct SVM
        # library(e1071)
        # svm.linear.res = svm(train.dat,train.response,type="C-classification",kernel="linear",cost=10)
        # svm.rbf.res = svm(train.dat,train.response,type="C-classification",kernel="radius",gamma=1,cost=10)
      }
    }
  }
  
}

