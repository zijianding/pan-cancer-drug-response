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
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/miRNAseq/")
pat.info = read.table("cancer.patient.drug.response.miRNAseq.20150519.txt",header=T,sep="\t",quote="")
cisplatin.dat = read.table("cisplatin.miRNAseq.gdac_20141206.preprocess.txt",header=T,sep="\t",quote="",row.names=1)


###data preprocess###
#patient information preprocess#
pat.info = pat.info[as.character(pat.info$data)=="YES.miRNAseq_mature.RPM_log2",]
#View(pat.info)
tmp.arr = vector(mode="character",length=nrow(pat.info))
tmp.arr[which(as.character(pat.info$response)=="Clinical Progressive Disease")] = "insensitive"
tmp.arr[which(pat.info$response=="Stable Disease")] = "insensitive"
tmp.arr[which(pat.info$response=="Complete Response")] = "sensitive"
tmp.arr[which(pat.info$response=="Partial Response")] = "sensitive"
pat.info$response = tmp.arr
cisplatin.info = pat.info[as.character(pat.info$drug)=="Cisplatin",]


#annotate patients in info
cancer_type = unique(as.character(cisplatin.info$cancer))
all_pats = vector(length=nrow(cisplatin.info),mode="character")
for(i in 1:length(all_pats))
{
  curr_pat = as.character(cisplatin.info$patient[i])
  curr_cancer = as.character(cisplatin.info$cancer[i])
  pat_arr = unlist(strsplit(curr_pat,split="\\-"))
  all_pats[i] = paste(curr_cancer,pat_arr[1],pat_arr[2],pat_arr[3],"01",sep=".")
}

sum(all_pats %in% colnames(cisplatin.dat)) < length(all_pats)

pat_ix = match(colnames(cisplatin.dat),all_pats)

all_pats = all_pats[pat_ix]
cisplatin.info = cisplatin.info[pat_ix,]




###shuffle response in each cancer and 
#partition every cancer type into 5-fold###
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/shuffle/shuffle_response.in_cancer/miRNAseq/")
for( shuffle.times in 1:100 )
{
  #shuffle response
  shuffle.info = cisplatin.info
  for(i in 1:length(cancer_type))
  {
    curr.info = shuffle.info[shuffle.info$cancer==cancer_type[i],]
    set.seed(i)
    shuffle.resp = sample(as.character(curr.info$response),replace=F)
    curr.info$response = shuffle.resp
    shuffle.info[shuffle.info$cancer==cancer_type[i],] = curr.info
  }
  
  #start to partition
  all_res = list()
  for( all in 1:20)
  {
    
    k=5
    partition_res = matrix("NULL",nrow=1,ncol=5+3)
    partition_pre = matrix( 0,ncol=1,nrow=nrow(shuffle.info) )
    
    #cancer type
    mytable = table(shuffle.info$cancer)
    mytable = mytable[mytable>0]
    cancer_type = names(mytable)
    
    #start to partition, meanwhile keep ratio of sensitive and insensitive
    for( i in 1:length(cancer_type) )
    {
      curr_info = shuffle.info[as.character(shuffle.info$cancer)==cancer_type[i],]
      pat_num = nrow(curr_info)
      sen_pats = as.character(curr_info$patient[curr_info$response=="sensitive"])
      insen_pats = as.character(curr_info$patient[curr_info$response=="insensitive"])
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
    #store current result
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
  
  ##output the shuffle results
  
  file_name = paste("cisplatin.miRNAseq_shuffle_info.",shuffle.times,".txt",sep="")
  write.table(partition_res,file_name,col.names=T,row.names=F,sep="\t",quote=F)
  
  
  
}



###etc file for qsub###
shuffle.times = 1:100
test_partition = 1:100
file_name = paste("/data/home/zding/drug_sensitivity/data/shuffle/shuffle_response.in_cancer/miRNAseq/cisplatin.miRNAseq_shuffle_info.",shuffle.times,".txt",sep="")


test = rep( test_partition,times=length(test_partition) )
shuffle = rep( shuffle.times,each=length(test_partition))
files = rep(file_name,each=length(test_partition))


df = data.frame(test=test,shuffle=shuffle,files=files)
write.table(df,"C:/Users/zding/workspace/projects/drug_sensitivity/data/shuffle/shuffle_response.in_cancer/miRNAseq/shuffle_etc.txt",
            col.names=F,row.names=F,quote=F,sep="\t")



