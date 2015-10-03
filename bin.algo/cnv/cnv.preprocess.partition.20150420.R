###functions###
partition_data <- function(obs, k){
  #k-fold partition of obs
  cv.mat = matrix("NULL",nrow=length(obs),ncol=k)
  
  fd.size = floor(length(obs)/k)
  tmp = length(obs)%%k
  if(tmp>0)
  {
    fd.ix = rep(1:tmp,each = (fd.size+1))
    fd.ix = c(fd.ix,rep((tmp+1):k,each=fd.size))
  }
  if(tmp==0)
  {
    fd.ix = rep(1:k,each=fd.size)
  }
  
  for( j in 1:k)
  {
    cv.mat[fd.ix!=j,j] = "train"
    cv.mat[fd.ix==j,j] = "validation"
  }
  rownames(cv.mat) = obs
  return(cv.mat)
}

partition_data_ass <- function(obs,info)
{
  ix = match(obs,info$patient)
  cancers = as.character( info$cancer[ix] )
  response = as.character( info$response[ix]  )
  mat = matrix("NULL",nrow=length(obs), ncol=3 )
  mat[,1] = obs
  mat[,2] = cancers
  mat[,3] = response
  rownames(mat) = mat[,1]
  return(mat)
}
###load data###
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv")
pat.info = read.table("cancer.patient.drug.response.gistic2.20150330.txt",header=T,sep="\t",quote="")
#cisplatin
#cisplatin.dat = read.table("Cisplatin.gistic2.gdac_20141206.txt",header=T,sep="\t",quote="",row.names=1)
#carboplatin
cisplatin.dat = read.table("Carboplatin.gistic2.gdac_20141206.txt",header=T,sep="\t",quote="",row.names=1)
#paclitaxel
#cisplatin.dat = read.table("Paclitaxel.gistic2.gdac_20141206.txt",header=T,sep="\t",quote="",row.names=1)

###data preprocess###
#patient information preprocess#
pat.info = pat.info[as.character(pat.info$data)=="YES.Gistic2_focal.continous",]
#View(pat.info)
tmp.arr = vector(mode="character",length=nrow(pat.info))
tmp.arr[which(as.character(pat.info$response)=="Clinical Progressive Disease")] = "insensitive"
tmp.arr[which(pat.info$response=="Stable Disease")] = "insensitive"
tmp.arr[which(pat.info$response=="Complete Response")] = "sensitive"
tmp.arr[which(pat.info$response=="Partial Response")] = "sensitive"
pat.info$response = tmp.arr
cisplatin.info = pat.info[as.character(pat.info$drug)=="Carboplatin",]
#data preprocess
#delete cancers with no more than 2 samples
mytable = table(cisplatin.info$cancer)
mytable = mytable[mytable<=2]
delete_cancer = names(mytable[mytable>0])
delete_ix = c()
for(i in 1:length(delete_cancer))
{
  ix = which(delete_cancer[i]==as.character(cisplatin.info$cancer))
  delete_ix = c(delete_ix, ix)
}
if(length(delete_ix)>0){
  cisplatin.info = cisplatin.info[-delete_ix,]
  
  ###find common patients###
  dat.pats = c()
  both.ix = c()
  for(i in 1:ncol(cisplatin.dat))
  {
    tmp = strsplit(as.character(colnames(cisplatin.dat)[i]),"\\.")
    tmp.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
    if( (!is.na( pmatch("01",tmp[[1]][5]) )) && (!is.na(match(tmp.pat,as.character(cisplatin.info$patient)))) )
    {
      #primary.tumor.arr = c(primary.tumor.arr,i)
      dat.pats = c(dat.pats,tmp.pat)
      both.ix = c(both.ix,i)
    }
    else
    {
      print(tmp[[1]][5])
    }
  }
  cisplatin.dat = cisplatin.dat[,both.ix]
  cis.info.delete.ix = c()
  for(i in 1:nrow(cisplatin.info))
  {
    tmp.pat = as.character(cisplatin.info$patient[i])
    if( is.na(match(tmp.pat,dat.pats)))
    {
      cis.info.delete.ix = c(cis.info.delete.ix,i)
    }
  }
  if(length(cis.info.delete.ix)>0)
  {
    cisplatin.info = cisplatin.info[-cis.info.delete.ix,]
  }
  
  
  #after deleting some samples, any newly "2 sample cancer"?
  mytable = table(cisplatin.info$cancer)
  mytable = mytable[mytable<=2]
  delete_cancer = names(mytable[mytable>0])
  delete_ix = c()
  for(i in 1:length(delete_cancer))
  {
    ix = which(delete_cancer[i]==as.character(cisplatin.info$cancer))
    delete_ix = c(delete_ix, ix)
  }
  if(length(delete_ix)>0){
    cisplatin.info = cisplatin.info[-delete_ix,]
    
    ###find common patients###
    dat.pats = c()
    both.ix = c()
    for(i in 1:ncol(cisplatin.dat))
    {
      tmp = strsplit(as.character(colnames(cisplatin.dat)[i]),"\\.")
      tmp.pat = paste(tmp[[1]][2],tmp[[1]][3],tmp[[1]][4],sep="-")
      if( (!is.na( pmatch("01",tmp[[1]][5]) )) && (!is.na(match(tmp.pat,as.character(cisplatin.info$patient)))) )
      {
        #primary.tumor.arr = c(primary.tumor.arr,i)
        dat.pats = c(dat.pats,tmp.pat)
        both.ix = c(both.ix,i)
      }
      else
      {
        print(tmp[[1]][5])
      }
    }
    cisplatin.dat = cisplatin.dat[,both.ix]
    cis.info.delete.ix = c()
    for(i in 1:nrow(cisplatin.info))
    {
      tmp.pat = as.character(cisplatin.info$patient[i])
      if( is.na(match(tmp.pat,dat.pats)))
      {
        cis.info.delete.ix = c(cis.info.delete.ix,i)
      }
    }
    if(length(cis.info.delete.ix)>0)
    {
      cisplatin.info = cisplatin.info[-cis.info.delete.ix,]
    }
    
  }
}





#delete genes that have zeros in more than 95% samples in each cancer
#NAs distribution
zero_num = c()
for( i in 1:nrow(cisplatin.dat))
{
  curr_na = sum( as.numeric(as.character(cisplatin.dat[i,])) == 0  )
  zero_num = c(zero_num,curr_na)
}
hist(zero_num/ncol(cisplatin.dat),100)
#
#ratio = 1 - quantile(zero_num/ncol(cisplatin.dat))[4]
ratio = 0.1
delete.ix = c()
for(i in 1:nrow(cisplatin.dat))
{
  tmp.dat = as.numeric(as.character(cisplatin.dat[i,]))
  tmp.zero = sum(tmp.dat==0)
  tmp.ratio = 1 - tmp.zero/ncol(cisplatin.dat)
  if( tmp.ratio<=ratio)
  {
    delete.ix = c(delete.ix,i)
  }
}
cisplatin.dat = cisplatin.dat[-delete.ix,]



###core cancer###
table.all_cancer = table(cisplatin.info$cancer)
df.all_cancer = data.frame(table.all_cancer)
df.all_cancer = df.all_cancer[-which(df.all_cancer$Freq==0),]
neg_ratio = vector(length=nrow(df.all_cancer),mode="numeric")
insen_num = vector(length=nrow(df.all_cancer),mode="numeric")
for(i in 1:nrow(df.all_cancer))
{
  df.tmp = cisplatin.info[as.character(cisplatin.info$cancer)==df.all_cancer[i,1],]
  neg_num = sum(df.tmp$response=="sensitive")
  neg_ratio[i] = neg_num/nrow(df.tmp)
  insen_num[i] = sum(as.character(df.tmp$response)=="insensitive")
}
coeff = df.all_cancer$Freq/neg_ratio
df.all_cancer = data.frame(df.all_cancer,neg_ratio,insen_num,coeff)
colnames(df.all_cancer) = c("cancer","observations","sensitive.ratio","insen_num","coeff")
df.all_cancer = df.all_cancer[order(df.all_cancer[,5],decreasing=T),]
View(df.all_cancer)


####sample test data#####

##partition every cancer type into 5-fold##
#repeat the process for 20 times##
all_res = list()
for(all in 1:20)
{
  k=5
  partition_res = matrix("NULL",nrow=1,ncol=5+3)
  partition_pre = matrix( 0,ncol=1,nrow=nrow(cisplatin.info) )
  
  #cancer type
  mytable = table(cisplatin.info$cancer)
  mytable = mytable[mytable>0]
  cancer_type = names(mytable)
  
  #start to partition, meanwhile keep ratio of sensitive and insensitive
  for( i in 1:length(cancer_type) )
  {
    curr_info = cisplatin.info[as.character(cisplatin.info$cancer)==cancer_type[i],]
    pat_num = nrow(curr_info)
    sen_pats = as.character(curr_info$patient[curr_info$response=="sensitive"])
    sen_pats = sample(sen_pats,size=length(sen_pats),replace=F)
    insen_pats = as.character(curr_info$patient[curr_info$response=="insensitive"])
    insen_pats = sample(insen_pats,size=length(insen_pats),replace=F)
    if( length(insen_pats)>=5 )
    {
      sen_mat = partition_data(sen_pats,k)
      insen_mat = partition_data(insen_pats,k)
      pats_mat = rbind(sen_mat,insen_mat)
    }
    if( length(insen_pats) < 5 )
    {
      pats = c(sen_pats,insen_pats)
      pats = sample(pats,replace=F)
      pats_mat = partition_data(pats,k)
      
    }
    #add information
    mat_info = partition_data_ass(rownames(pats_mat),curr_info)
    mat_res = cbind(mat_info,pats_mat) 
    partition_res = rbind(partition_res,mat_res)
  }
  partition_res = partition_res[-1,]
  for(i in 1:nrow(partition_res))
  {
    curr_pat = partition_res[i,1]
    curr_cancer = partition_res[i,2]
    tmp = unlist(strsplit(curr_pat,split="\\-"))
    tmp = paste(curr_cancer,tmp[1],tmp[2],tmp[3],"01",sep=".")
    partition_res[i,1] = tmp
  }
  colnames(partition_res) = c("patient","cancer","response",as.character(seq(1,k,by=1)))
  
  all_res[[all]] = partition_res
  
}

#merge all the partition results
#next sentences may be wrong
partition_res = matrix("NULL",nrow=1,ncol=100+3)
pats = all_res[[1]][,1]
rownames(all_res[[1]]) = NULL
for(all in 2:20)
{
  all_res[[all]] = all_res[[all]][match(pats,all_res[[all]][,1]),]
  all_res[[all]] = all_res[[all]][,-c(1,2,3)]
  rownames(all_res[[all]]) = NULL
}
partition_res = do.call(cbind,all_res)





###output###
#deal with colnames
tmp = colnames(cisplatin.dat)
colnames(cisplatin.dat) = substring(tmp,1,nchar(tmp)-1)
write.table(cisplatin.dat,"cisplatin.gistic2_focal.gdac_20141206.preprocess.txt",quote=F,row.names=T,col.names=T,sep="\t")

write.table(partition_res,"cisplatin.gistic2_focal.5_fold_cv.mat.txt",col.names=T,row.names=F,sep="\t",quote=F)




