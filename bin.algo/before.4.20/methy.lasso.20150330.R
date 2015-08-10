#setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/methylation")
args <- commandArgs(trailingOnly=TRUE)
###patient information
#patient.info =  read.table("cancer.patient.drug.response.methylation.20150303.txt",header=T,sep="\t",quote="")
patient.info = read.table(args[1],header=T,sep="\t",quote="")
patient.info$response = as.character(patient.info$response)

patient.info = patient.info[as.character(patient.info$data)=="YES.methylation",]
patient.info$response[as.character(patient.info$response)=="Clinical Progressive Disease"]="insensitive"
patient.info$response[patient.info$response=="Stable Disease"] = "insensitive"
patient.info$response[patient.info$response=="Complete Response"] = "sensitive"
patient.info$response[patient.info$response=="Partial Response"] = "sensitive"



#############################Cisplatin#############################
cisplatin.info = patient.info[as.character(patient.info$drug)=="Cisplatin",]

######preprocess methylation data################
#patient.data = read.table("Cisplatin.methylation450k.gdac_20141206.preprocess.txt",sep="\t",header=T,quote="",row.names=1)
patient.data = read.table(args[2],sep="\t",header=T,quote="",row.names=1)
primary.tumor.arr = c()
for( i in 2:ncol(patient.data) )
{
  tmp = strsplit(as.character(colnames(patient.data)[i]),"\\.")
  cancer = tmp[[1]][1]
  curr.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
  if( tmp[[1]][5] == "01A" )
  {
    primary.tumor.arr = c(primary.tumor.arr,i)
  }
  else
  {
    print(tmp[[1]][5])
  }
}

patient.data = patient.data[,primary.tumor.arr] #then no column of gene
for(i in 1:ncol(patient.data))
{
  patient.data[,i] = as.numeric(as.character(patient.data[,i]))
}

#patients divided to "sensitive" and "insensitive" groups

sen.index = c()
insen.index = c()
patients = colnames(patient.data)
for(i in 1:length(patients))
{
  tmp = strsplit(patients[i],"\\.")
  cancer = tmp[[1]][1]
  curr.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
  pat.ix = match(curr.pat,as.character(cisplatin.info$patient))
  
  if( cisplatin.info$response[pat.ix] == "sensitive")
  {
    sen.index = c(sen.index,i )
  }
  if( cisplatin.info$response[pat.ix] == "insensitive")
  {
    insen.index = c(insen.index,i)
  }
}


####start to find biomarkers#####


#t test
test.res = matrix(0,nrow=nrow(patient.data),ncol=4)
colnames(test.res) = c("sensitive","insensitive","pvalue","fdr")
rownames(test.res) = rownames(patient.data)


for( i in 1:nrow(patient.data) )
{
  sen.data = unlist(patient.data[i,sen.index])
  insen.data = unlist(patient.data[i,insen.index])
  curr.res = t.test(sen.data,insen.data)
  test.res[i,3] = curr.res$p.value
  test.res[i,1] = mean(sen.data)
  test.res[i,2] = mean(insen.data)
}

test.res[,4] = p.adjust(test.res[,3],method="fdr")
sum(test.res[,4]<=0.05)
#test.res = test.res[order(test.res[,4],decreasing=F),]
#patient.data = patient.data[order(test.res[,4],decreasing=F),]

#select the probes, also in the train data
pre.probe.ix = which(test.res[,4]<=0.05)


cisplatin.data = patient.data[pre.probe.ix,]





###start to calculate the single cancer difference
cancers = unique(as.character(cisplatin.info$cancer))
#return the cancer index in the data
cancer.index = vector(mode="character",length=ncol(patient.data))
patients = as.character(colnames(patient.data))
for( i in 1:ncol(patient.data))
{
  tmp = strsplit(patients[i],"\\.")
  cancer.index[i] = tmp[[1]][1]
}

res.mat = matrix(0,nrow(cisplatin.data),ncol=2*length(cancers))
for( i in 1:nrow(cisplatin.data))
{
  for( j in 1:length(cancers))
  {
    curr.index = which(cancers[j]==cancer.index)
    curr.sen.ix = intersect(curr.index,sen.index)
    curr.insen.ix = intersect(curr.index,insen.index)
    sen.data = unlist(cisplatin.data[i,curr.sen.ix])
    insen.data = unlist(cisplatin.data[i,curr.insen.ix])
    
    if((length(sen.data)==0 ) && (length(insen.data)>=1) )
    {
      res.mat[i,2*j-1] =  - mean(insen.data)
      res.mat[i,2*j] = NA
    }
    if( (length(sen.data)>=1 ) && (length(insen.data)==0) )
    {
      res.mat[i,2*j-1] =  mean(sen.data)
      res.mat[i,2*j] = NA
    }
    if( (length(sen.data)>=1 ) && (length(insen.data)>=1) )
    {
      res.mat[i,2*j-1] = mean(sen.data) - mean(insen.data)
      res.mat[i,2*j] = NA
    }
    if( (length(sen.data)>1) &&  (length(insen.data)>1) ) 
    {
      curr.res = t.test(sen.data,insen.data)
      res.mat[i,2*j-1] = mean(sen.data) - mean(insen.data)
      res.mat[i,2*j] = curr.res$p.value
    } 
  }
}
rownames(res.mat) = rownames(cisplatin.data)
colnames(res.mat) = rep(c("sen/insen","t.test.p.value"),14)


response = vector(mode="character",length=ncol(patient.data))
response[sen.index] = "sensitive"
response[insen.index] = "insensitive"



#arrange data
#probe -> Label(sen/insen) -> variable(tumor/pancancer) -> value
pdf("probe.val_distribution.pancancer.20150331.pdf")
p.value.order = order(test.res[pre.probe.ix,4])
require(ggplot2)
for(i in 1:length(p.value.order) )
{
  #for current probe
  
  curr.dat = patient.data[pre.probe.ix[p.value.order[i]],]
  probe = rownames(curr.dat)
  curr.dat = as.numeric(curr.dat[1,])
  curr.dat = c(curr.dat,curr.dat)
  label=as.vector(response)
  label = c(label,label)
  variable = c( as.vector(cancer.index),rep("ALL",length(cancer.index)) )
  df = data.frame(label,variable,curr.dat)
  colnames(df) = c("Label","variable","beta.value")
  ggplot(data=df,aes(x=as.factor(variable),y=beta.value)) + geom_boxplot(aes(fill=Label),notch=F) + labs(title=paste(probe,"across cancer types",sep=" "),x="cancer type",y="beta value")
}
dev.off()































#estimate NAs using mean values
for(i in 1:ncol(cisplatin.data))
{
  tmp = which(is.na(cisplatin.data[,i]))
  if( length(tmp)>0 )
  {
    curr.mean = mean(cisplatin.data[,i],na.rm=T)
    cisplatin.data[tmp,i] = curr.mean
  }
}

#data


response.train = as.factor(response[train.ix])
response.test = response[valid.ix]


#start to fit model, find the lambda
library(glmnet)
#fit.lasso = glmnet(x=cisplatin.data,y=response.train,family="binomial",alpha=1)
cv.lasso = cv.glmnet(x=cisplatin.data,y=response.train,alpha=1,nfolds=5,type.measure="class",family="binomial")

plot(cv.lasso)
best.lam = cv.lasso$lambda.min
cisplatin.test = t(patient.data[pre.probe.ix,valid.ix])
for(i in 1:ncol(cisplatin.test))
{
  tmp = which(is.na(cisplatin.test[,i]))
  if( length(tmp)>0 )
  {
    curr.mean = mean(cisplatin.test[,i],na.rm=T)
    cisplatin.test[tmp,i] = curr.mean
  }
}
glm.lasso = glmnet(x=cisplatin.data,y=response.train,lambda=cv.lasso$lambda.min,family="binomial",alpha=1)
prediction.res = predict(object = glm.lasso,newx = cisplatin.test,s=cv.lasso$lambda.min,type="response")
prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
prediction.res.arr[which(prediction.res>=0.5)] = levels(response.train)[2]
prediction.res.arr[which(prediction.res<0.5)] = levels(response.train)[1]

mis.rate = sum(prediction.res.arr!=response.test)
mis.rate = mis.rate/length(prediction.res.arr)
mis.rate
prediction.res.arr[which(prediction.res.arr!=response.test)]

####several question#####
#misclassified evenly distributed across cancer types?
#correctly classified evenly distributed across cancer types?
  






##################single cancer###############
cancers = unique(cisplatin.info$cancer)
#TGCT as train, LUAD as test
patients = as.character(colnames(patient.data))
tgct.ix = c()
luad.ix = c()
blca.ix = c()
for( i in 1:length(patients))
{
  tmp = strsplit(patients[i],"\\.")
  cancer = tmp[[1]][1]
  
  if( cancer == "TGCT" )
  {
    tgct.ix = c(tgct.ix,i)
  }
  if( cancer == "LUAD")
  {
    luad.ix = c(luad.ix,i)
  }
  if( cancer == "BLCA")
  {
    blca.ix = c(blca.ix,i)
  }
}

tgct.dat = patient.data[,tgct.ix]
luad.dat = patient.data[,luad.ix]
blca.dat = patient.data[,blca.ix]

response.tgct = response[tgct.ix]#1:40/insenstive:sensitive
response.luad = response[luad.ix]#8:32
response.blca = response[blca.ix]#14:20


#train data: luad; test data: blca
luad.test = matrix(0,nrow=nrow(luad.dat),ncol=4)
colnames(luad.test) = c("sensitive.mean","insensitive.mean","p.value","bonferroni")
rownames(luad.test) = rownames(luad.dat)
for(i in 1:nrow(luad.dat))
{
  sen.dat = luad.dat[,response.luad=="sensitive"]
  insen.dat = luad.dat[,response.luad=="insensitive"]
  curr.res = t.test(sen.data,insen.data)
  test.res[i,3] = curr.res$p.value
  test.res[i,1] = mean(sen.data)
  test.res[i,2] = mean(insen.data)
  
}
test.res[,4] = p.adjust(test.res[,3],method="bonferroni")
print(sum(test.res[,4]<=0.05))

library(glmnet)
pre.probe.ix = which(test.res[,4]<=0.05)
luad.dat1 = t(luad.dat[pre.probe.ix,])

#estimate NAs using mean values
for(i in 1:ncol(luad.dat1))
{
  tmp = which(is.na(luad.dat1[,i]))
  if( length(tmp)>0 )
  {
    curr.mean = mean(luad.dat1[,i],na.rm=T)
    luad.dat1[tmp,i] = curr.mean
  }
}

cv.lasso = cv.glmnet(x=luad.dat1,y=as.factor(response.luad),alpha=1,nfolds=5,type.measure="class",family="binomial")
#plot(cv.lasso)
best.lam = cv.lasso$lambda.min
glm.lasso = glmnet(x=luad.dat1,y=as.factor(response.luad),lambda=cv.lasso$lambda.min,family="binomial",alpha=1)

###test on BLCA
blca.dat1 = t(blca.dat1[pre.probe.ix,])
#estimate NAs using mean values
for(i in 1:ncol(blca.dat1))
{
  tmp = which(is.na(blca.dat1[,i]))
  if( length(tmp)>0 )
  {
    curr.mean = mean(blca.dat1[,i],na.rm=T)
    blca.dat1[tmp,i] = curr.mean
  }
}

prediction.res = predict(object = glm.lasso,newx = blca.dat1,s=cv.lasso$lambda.min,type="response")
prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(response.luad))[2]
prediction.res.arr[which(prediction.res<0.5)] = levels(response.luad)[1]

mis.rate = sum(prediction.res.arr!=response.blca)
mis.rate = mis.rate/length(prediction.res.arr)
mis.rate
#prediction.res.arr[which(prediction.res.arr!=response.test)]

###test on tgct
tgct.dat1 = t(tgct.dat1[pre.probe.ix,])
#estimate NAs using mean values
for(i in 1:ncol(tgct.dat1))
{
  tmp = which(is.na(tgct.dat1[,i]))
  if( length(tmp)>0 )
  {
    curr.mean = mean(tgct.dat1[,i],na.rm=T)
    tgct.dat1[tmp,i] = curr.mean
  }
}

prediction.res = predict(object = glm.lasso,newx = tgct.dat1,s=cv.lasso$lambda.min,type="response")
prediction.res.arr = vector(mode="character",length=nrow(prediction.res))
prediction.res.arr[which(prediction.res>=0.5)] = levels(as.factor(response.luad))[2]
prediction.res.arr[which(prediction.res<0.5)] = levels(response.luad)[1]

mis.rate = sum(prediction.res.arr!=response.tgct)
mis.rate = mis.rate/length(prediction.res.arr)
mis.rate
#prediction.res.arr[which(prediction.res.arr!=response.test)]



