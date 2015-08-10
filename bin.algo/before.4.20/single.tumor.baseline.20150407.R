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
##################
#draw the cancer.df
# tmp.cancer.arr = rep(rownames(cancer.df),each=2)
# tmp.resp.arr = rep(c("sensitive","insensitive"),nrow(cancer.df))
# tmp.arr = vector(length=2*nrow(cancer.df),mode="numeric")
# tmp.arr[seq(1,(length(tmp.arr)-1),by=2)] = cancer.df$sensitive
# tmp.arr[seq(2,length(tmp.arr),by=2)] = cancer.df$insensitive
# res.df2 = data.frame(tmp.cancer.arr,tmp.resp.arr,tmp.arr)
# colnames(res.df2) = c("cancer","response","patients")
# library(ggplot2)
# library(plyr)
# res.df2 <- ddply(res.df2, .(cancer), transform, pos = cumsum(patients) - (0.5 * patients))
# 
# ggplot(res.df2, aes(x = cancer, y = patients)) +
#   geom_bar(aes(fill = response), stat="identity") +
#   geom_text(aes(label = patients, y = pos), size = 3) +
#   labs(title="Clinical drug responses of patients across cancer types",x="cancer type",y="Patients") +
#   theme(text=element_text(size=16),plot.title=element_text(size=20,face="bold"))
cancer.df = cancer.df[cancer.df[,1]!=1,]
#View(cancer.df)
#library(ggplot2)
#ggplot(cancer.df,aes(x=colnames(cancer.df),y=patient)) + geom_bar(aes(fill=sensitive),stat="identity")
#core cancer: LUAD, BLCA, MESO, STAD

#########################
####model constrution####
core.cancer = c("LUAD","BLCA","MESO","STAD")
library(glmnet)
#library(e1071)
for( r in 1:length(core.cancer) )
{
  tmp.char = paste("The current tumor is",core.cancer[r],sep=" ")
  print(tmp.char)
  #############################
  #find samples of current tumor
  patients = colnames(patient.data)
  curr.sen = c()
  curr.insen = c()
  pat.ix = c()
  for(i in 1:ncol(patient.data) )
  {
    tmp = strsplit(patients[i],"\\.")
    cancer = tmp[[1]][1]
    curr.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
    curr.ix = match(curr.pat,as.character(cisplatin.info$patient))
    if( cancer == core.cancer[r] )
    {
      #pat.ix = c(pat.ix,curr.ix)
      pat.ix = c(pat.ix,i) #the index of tumor data in all data
      response = as.character(cisplatin.info$response[curr.ix])
      if( response == "sensitive" )
      {
        curr.sen = c(curr.sen,length(pat.ix))
      }
      if( response == "insensitive" )
      {
        curr.insen = c( curr.insen, length(pat.ix) )
      }
    }
  }
  curr.dat = patient.data[,pat.ix]
  response.arr = vector(length=length(pat.ix),mode="character")
  response.arr[curr.sen] = "sensitive"
  response.arr[curr.insen] = "insensitive"
  #########################
  #5-fold cross validation
  #fd.ix = vector(length=ncol(curr.dat),mode="numeric")
  fd.size = floor(ncol(curr.dat)/5)
  tmp = ncol(curr.dat)%%5
  if(tmp>0)
  {
    fd.ix = rep(1:tmp,each = (fd.size+1))
    fd.ix = c(fd.ix,rep((tmp+1):5,each=fd.size))
  }
  if(tmp==0)
  {
    fd.ix = rep(1:5,each=fd.size)
  }
  
  #fd.ix = c(fd.ix,rep(5,tmp))
  for( f in 1:5 )
  {
    tmp.char = paste("The current test fold of",core.cancer[r],"is",f,sep=" ")
    print(tmp.char)
    ############################
    train.ix = which(fd.ix != f)
    test.ix = which(fd.ix == f )
    #######
    #t test
    pval.arr = vector(mode="numeric",length=nrow(curr.dat))
    for(j in 1:nrow(curr.dat))
    {
      sen.dat = unlist(curr.dat[j,intersect(train.ix,curr.sen)])
      insen.dat = unlist(curr.dat[j,intersect(train.ix,curr.insen)])
      test.res = t.test(sen.dat,insen.dat)
      #pval.arr = c(pval.arr,test.res$p.value)
      pval.arr[j] = test.res$p.value
    }
    fdr.arr = p.adjust(pval.arr,method="fdr")
    ##################################
    #if there is no significant probes
    if( sum(fdr.arr<=0.05)<=1 )
    {
      tmp.char = paste("NO sig probes in fold",f,"of tumor",core.cancer[r],sep=" ")
      print(tmp.char)
      #the 1st 1000:10000,pace=10000 probes
      fdr.arr.ix = sort(fdr.arr,index.return=T)$ix
      probe.num = seq(1000,10000,by=1000)
      for( p in 1:length(probe.num) )
      {
        tmp.char = paste("The first",probe.num[p],"probes","in fold",f,"of tumor",core.cancer[r],sep=" ")
        print(tmp.char)
        ########################################
        pre.probe.ix = fdr.arr.ix[1:probe.num[p]]
        #now.dat = curr.dat[pre.probe.ix,]
        ###############################
        #training data and testing data
        train.dat = t(curr.dat[pre.probe.ix,train.ix])
        for(i in 1:ncol(train.dat))
        {
          tmp = which(is.na(train.dat[,i]))
          if( length(tmp)>0 )
          {
            curr.mean = mean(train.dat[,i],na.rm=T)
            train.dat[tmp,i] = curr.mean
          }
        }
        train.response = response.arr[train.ix]
        test.dat = t(curr.dat[pre.probe.ix,test.ix])
        for(i in 1:ncol(test.dat))
        {
          tmp = which(is.na(test.dat[,i]))
          if( length(tmp)>0 )
          {
            curr.mean = mean(test.dat[,i],na.rm=T)
            test.dat[tmp,i] = curr.mean
          }
        }
        test.response = response.arr[test.ix]
        ##############################
        #construct logistic regression
        cv.res = cv.glmnet(x=train.dat,y=as.factor(train.response),alpha=1,nfolds=5,type.measure="class",family="binomial")
        plot(cv.res)
        
        glm.res = glmnet(x=train.dat,y=as.factor(train.response),lambda=cv.res$lambda.min,family="binomial",alpha=1)
        
        prediction.res = predict(object = glm.res,newx = test.dat,s=cv.res$lambda.min,type="response")
        prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
        prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(train.response))[2]
        prediction.res.arr[which(prediction.res<0.5)] = levels(as.factor(train.response))[1]
        
        
        #mis.rate = sum(prediction.res.arr!=test.response)
        #mis.rate = mis.rate/length(prediction.res.arr)
        ###the proprotion
        insen.ratio = sum(test.response=="insensitive")/length(test.response)
        tmp.char = paste("The proportion of insensitive samples is", insen.ratio,sep=" ")
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
          tmp.char = paste("The FDR(small)",fdr,"of fold",f,"in",core.cancer[r],sep=" ")
          print(tmp.char)
        }
        if((tp+fp)==0)
        {
          tmp.char = paste("The FDR(small) is NA of fold",f,"in",core.cancer[r],sep=" ")
          print(tmp.char)
        }
        #Recall = TP/(TP+FN), as large as possible
        if((tp+fn)>0)
        {
          tmp.char = paste("The recall(large) is",tp/(tp+fn),"of fold",f,"in",core.cancer[r],sep=" ")
          print(tmp.char)
        }
        if((tp+fn)==0)
        {
          tmp.char = paste("The recall(large) is NA of fold",f,"in",core.cancer[r],sep=" ")
          print(tmp.char)
        }
        #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
        if((fp+tn)>0)
        {
          tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"of fold",f,"in",core.cancer[r],sep=" ")
          print(tmp.char)
        }
        if((fp+tn)==0)
        {
          tmp.char = paste("The FPR(small) is NA of fold",f,"in",core.cancer[r],sep=" ")
          print(tmp.char)
        }
        ###the misclassification rate
        mis.rate = (fp+fn)/(fp+fn+tn+tp)
        tmp.char = paste("The misclassification rate is",mis.rate,sep=" ")
        print(tmp.char)
        
        #construct SVM
#         library(e1071)
#         svm.linear.res = svm(train.dat,train.response,type="C-classification",kernel="linear",cost=10)
#         svm.rbf.res = svm(train.dat,train.response,type="C-classification",kernel="radius",gamma=1,cost=10)
        
        
      }
      
    }
    #if there is any significant probes
    else
    {
      sig.num = sum(fdr.arr<=0.05)
      tmp.char = paste(sig.num,"sig probes of fold ",f,"in",core.cancer[r],sep=" ")
      print(tmp.char)
      
      pre.probe.ix = which(fdr.arr<=0.05)
      #curr.dat = curr.dat[pre.probe.ix,]
      
      train.dat = t(curr.dat[pre.probe.ix,train.ix])
      for(i in 1:ncol(train.dat))
      {
        tmp = which(is.na(train.dat[,i]))
        if( length(tmp)>0 )
        {
          curr.mean = mean(train.dat[,i],na.rm=T)
          train.dat[tmp,i] = curr.mean
        }
      }
      train.response = response.arr[train.ix]
      test.dat = t(curr.dat[pre.probe.ix,test.ix])
      for(i in 1:ncol(test.dat))
      {
        tmp = which(is.na(test.dat[,i]))
        if( length(tmp)>0 )
        {
          curr.mean = mean(test.dat[,i],na.rm=T)
          test.dat[tmp,i] = curr.mean
        }
      }
      test.response = response.arr[test.ix]
      
      
      #construct logistic regression
      cv.res = cv.glmnet(x=train.dat,y=as.factor(train.response),alpha=1,nfolds=5,type.measure="class",family="binomial")
      plot(cv.res)
      
      
      glm.res = glmnet(x=train.dat,y=as.factor(train.response),lambda=cv.res$lambda.min,family="binomial",alpha=1)
      
      prediction.res = predict(object = glm.res,newx = test.dat,s=cv.res$lambda.min,type="response")
      prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
      prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(train.response))[2]
      prediction.res.arr[which(prediction.res<0.5)] = levels(as.factor(train.response))[1]
      
      
      
      #mis.rate = sum(prediction.res.arr!=test.response)
      #mis.rate = mis.rate/length(prediction.res.arr)
      ###the proprotion
      insen.ratio = sum(test.response=="insensitive")/length(test.response)
      tmp.char = paste("The proportion of insensitive samples is", insen.ratio,sep=" ")
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
        tmp.char = paste("The FDR(small)",fdr,"of fold",f,"in",core.cancer[r],sep=" ")
        print(tmp.char)
      }
      if((tp+fp)==0)
      {
        tmp.char = paste("The FDR(small) is NA of fold",f,"in",core.cancer[r],sep=" ")
        print(tmp.char)
      }
      #Recall = TP/(TP+FN), as large as possible
      if((tp+fn)>0)
      {
        tmp.char = paste("The recall(large) is",tp/(tp+fn),"of fold",f,"in",core.cancer[r],sep=" ")
        print(tmp.char)
      }
      if((tp+fn)==0)
      {
        tmp.char = paste("The recall(large) is NA of fold",f,"in",core.cancer[r],sep=" ")
        print(tmp.char)
      }
      #FPR=FP/(FP+TN),as small as possible, but may be TN is really large
      if((fp+tn)>0)
      {
        tmp.char = paste("The FPR(small) is ",fp/(fp+tn),"of fold",f,"in",core.cancer[r],sep=" ")
        print(tmp.char)
      }
      if((fp+tn)==0)
      {
        tmp.char = paste("The FPR(small) is NA of fold",f,"in",core.cancer[r],sep=" ")
        print(tmp.char)
      }
      ###the misclassification rate
      mis.rate = (fp+fn)/(fp+fn+tn+tp)
      tmp.char = paste("The misclassification rate is",mis.rate,sep=" ")
      print(tmp.char)
      
      #construct SVM
      #svm.linear.res = svm(train.dat,train.response,type="C-classification",kernel="linear",cost=10)
      #svm.rbf.res = svm(train.dat,train.response,type="C-classification",kernel="radius",gamma=1,cost=10)
    }
    
  }
  
}










