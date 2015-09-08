###load data###
args <- commandArgs(trailingOnly=TRUE)
#input
info_file = args[1]
input_folder = args[2] #no "/" at last
param_file_pattern = "pan.elanet.param.test_[0-9]*.20150701.txt"
roc_file_pattern = "pan.elanet.mid_res.test_[0-9]*.20150701.txt"
#output
output_folder = args[3]
create_folder = args[4]


source("source_all.R")

cisplatin.info = read.table(info_file,header=T,sep="\t",quote="")

dirs = list.dirs(path=input_folder,full.names=T,recursive=F)
output_names = lapply(dirs,
                      function(x){ tmp = strsplit(x,split="\\/"); 
                                   tmp=unlist(tmp); 
                                   return(tmp[length(tmp)]) } 
                      )
dir.create(file.path(output_folder,create_folder),showWarnings=FALSE)
setwd( file.path(output_folder,create_folder) )
for( j in 1:length(dirs))
{
  curr_folder = dirs[j]
  output_file_name = output_names[[j]]
  
  pan.files = list.files(path=curr_folder,full.names=T,pattern=param_file_pattern)
  pan.files = sort(pan.files)
  param = lapply(pan.files,read_table)
  filenames = lapply(pan.files,find_name)
  names( param ) = filenames
  
  
  score.files = list.files(path=curr_folder,full.names=T,pattern=roc_file_pattern)
  score.files = sort(score.files)
  pan_score = lapply(score.files,read_table)
  filenames = lapply(score.files,find_name)
  names(pan_score) = filenames
  
  #debug#
  if( length(param)!=length(pan_score))
  {
    cat("length of param is ",length(param),
        " and length of scores is ",length(pan_score),"\n")
    cat("The folder is",output_file_name,"\n")
    #next
  }
  
  ##calculate roc&auc##
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  roc_res <- foreach(i=1:length(pan_score) ) %dopar%
  {
    #for glmnet, second class is target class
    curr_pats = colnames(pan_score[[i]])
    resp = as.character( cisplatin.info$response[ match( curr_pats, as.character(cisplatin.info$patient) ) ] )
    #roc_res[[i]] = ensemble_roc(pan_score[[i]],resp,"sensitive")
    return( ensemble_roc(pan_score[[i]],resp,"sensitive") )
  }
  stopImplicitCluster()
  stopCluster(cl)
  auc_res = lapply( roc_res,function(x)return(trapz(x[,2],x[,1])) )
  auc_res = unlist(auc_res)
  
  
  ##get paramters##
  p.thresh = unlist( lapply(param,function(x)return(x[1,2])) )
  sig.gene = unlist( lapply(param,function(x)return(x[2,2])) )
  freq.thresh = unlist( lapply(param,function(x)return(x[3,2])) )
  freq.gene = unlist( lapply(param,function(x)return(x[4,2])) )
  
  
#   df = data.frame(p.thresh=p.thresh,sig.gene=sig.gene,
#                   freq.thresh=freq.thresh,freq.gene = freq.gene,
#                   auc=auc_res)
  
  #write.table(df,output_file_name,row.names=F,col.names=T,sep="\t",quote=F)
  #for debug
  write.table(auc_res,paste(output_file_name,"auc","txt",sep="."),
              row.names=F,col.names=T,sep="\n",quote=F)
}


