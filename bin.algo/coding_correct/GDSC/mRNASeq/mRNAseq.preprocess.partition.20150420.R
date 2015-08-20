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
#desktop
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/data/")
pat.info = read.table("cisplatin.ic50.txt",header=T,sep="\t",quote="")
cisplatin.dat = read.table("en_input_w5.txt",header=T,sep="\t",quote="",row.names=1)

###data preprocess###
##calibration data and patient information
dat.pats = c()
both.ix = c()
for(i in 1:ncol(cisplatin.dat))
{
  curr.pat = colnames(cisplatin.dat)[i]
  
  if( !is.na( match(curr.pat,pat.info$patient) )  )
  {
    #primary.tumor.arr = c(primary.tumor.arr,i)
    dat.pats = c(dat.pats,curr.pat)
    both.ix = c(both.ix,i)
  }
  else
  {
    print(curr.pat)
  }
}
cisplatin.dat = cisplatin.dat[,both.ix]
cis.info.delete.ix = c()
for(i in 1:nrow(pat.info))
{
  curr.pat = as.character(pat.info$patient[i])
  if( is.na(match(curr.pat,dat.pats)))
  {
    cis.info.delete.ix = c(cis.info.delete.ix,i)
  }
}
if(length(cis.info.delete.ix)>0)
{
  pat.info = pat.info[-cis.info.delete.ix,]
}
cisplatin.info = pat.info


#delete the tissue type
cisplatin.dat = cisplatin.dat[-((nrow(cisplatin.dat)-12):nrow(cisplatin.dat)),]

####data partition for each CV in 5-fold CV#####
###core cancer###
#define core cancer
tmp.df = data.frame(table(pat.info$cancer))
tmp.df = tmp.df[order(tmp.df[,2],decreasing=T),]
View(tmp.df)

####partition data into 10-fold#####
##partition every cancer type into 5-fold##
all_res = list()
for( all in 1:20)
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
    pats = as.character(curr_info$patient)

    #sen_mat = partition_data(sen_pats,k)
    #insen_mat = partition_data(insen_pats,k)
    #pats_mat = rbind(sen_mat,insen_mat)
    pats_mat = partition_data(pats,k)  
      
    #add information
    mat_info = partition_data_ass(rownames(pats_mat),curr_info)
    mat_res = cbind(mat_info,pats_mat) 
    partition_res = rbind(partition_res,mat_res)
  }
  partition_res = partition_res[-1,]
#   for(i in 1:nrow(partition_res))
#   {
#     curr_pat = partition_res[i,1]
#     curr_cancer = partition_res[i,2]
#     #tmp = unlist(strsplit(curr_pat,split="\\-"))
#     #tmp = paste(curr_cancer,tmp[1],tmp[2],tmp[3],"01",sep=".")
#     partition_res[i,1] = paste(curr_cancer,curr_pat)
#   }
  colnames(partition_res) = c("patient","cancer","response",as.character(seq(1,k,by=1)))
  all_res[[all]] = partition_res
}

#merge all the partition results
partition_res = matrix("NULL",nrow=nrow(all_res[[1]]),ncol=100+3)
pats = all_res[[1]][,1]
rownames(all_res[[1]]) = NULL
for(all in 2:20)
{
  all_res[[all]] = all_res[[all]][match(pats,all_res[[all]][,1]),]
  all_res[[all]] = all_res[[all]][,-c(1,2,3)]
  rownames(all_res[[all]]) = NULL
}
partition_res = do.call(cbind,all_res)


###output
write.table(cisplatin.dat,"cisplatin.mRNAseq.gdsc.preprocess.txt",quote=F,row.names=T,col.names=T,sep="\t")

write.table(partition_res,"cisplatin.mRNAseq_fold_cv.gdsc.txt",col.names=T,row.names=F,sep="\t",quote=F)

