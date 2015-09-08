###load data###
#source in
args <- commandArgs(trailingOnly=TRUE)
info_file = args[1]
real_folder = args[2]
shuffle_folder = args[3]
file_pattern = "pan.elanet.mid_res.test_[0-9]*.20150701.txt"
#source output
output_folder = args[4]
create_folder = args[5]
#functions
source("source_all.R")
#read
info = read.table(args[1],header=T,sep="\t",quote="")
real_files = list.files(path=read_folder,full.names=T,recursive=F,pattern=file_pattern)
real_scores = lapply(real_files,read_table)

shuffle_folders = list.dirs(path=shuffle_folder,full.names=T,recursive=F)

#calculate ROC and AUC#
#each curve
cl = makeCluster(no_cores)
registerDoParallel(cl)
real_roc <- foreach(i=1:length(real_scores)) %dopar%{
  curr_pats = colnames(real_scores[[i]])
  resp = as.character( info$response[ match( curr_pats, as.character(info$patient) ) ] )
  return( ensemble_roc(real_scores[[i]],resp,"sensitive") )
}
stopImplicitCluster()
stopCluster(cl)
#ave_curve
real_ave_roc = ave_roc(real_roc,type="vertical",title_str="real",measure="sd")

#shuffle curves
shuffle_ave = list()
for(j in 1:length(shuffle_folders))
{
  curr_folder = shuffle_folders[j]
  curr_files = list.files(path=curr_folder, full.names=T,recursive=F,pattern=file_pattern)
  curr_scores = lapply(curr_files,read_table)
  
  #all roc
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  curr_roc <- foreach(i=1:length(curr_scores)) %dopar%{
    curr_pats = colnames(curr_scores[[i]])
    resp = as.character( info$response[ match( curr_pats, as.character(info$patient) ) ] )
    return( ensemble_roc(curr_scores[[i]],resp,"sensitive") )
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  #ave roc
  shuffle_ave[[j]] = ave_roc(real_roc,type="vertial",title_str="shuffle",measure="sd")
}


#draw the distribution
shuffle_auc = sapply(shuffle_ave,function(x)return(x[[3]]))
df = data.frame(x=shuffl_auc)
real_auc = real_ave_roc[[3]]

fig = ggplot(data=df,aes(x=x)) + geom_histogram(alpha=0.3,colour='black',binwidth=0.1)
fig + geom_vline(xintercept=auc)
