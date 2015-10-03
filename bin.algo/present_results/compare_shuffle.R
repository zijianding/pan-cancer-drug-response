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
info = read.table(info_file,header=T,sep="\t",quote="")
real_files = list.files(path=real_folder,full.names=T,recursive=F,pattern=file_pattern)
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
real_auc = real_ave_roc[[3]][1]
cat("0\t",real_auc,"\n")

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
  curr_auc = ave_roc(curr_roc,type="vertical",title_str="shuffle",measure="sd")
  shuffle_ave[[j]] = curr_auc
  tmp1 = length(curr_auc)
  tmp2 = curr_auc[[3]][1]
  cat(j,"\t",tmp1,"\t",tmp2,"\n")
}


#draw the distribution


shuffle_auc = lapply(shuffle_ave,function(x){return(x[[3]][1])})
shuffle_auc = unlist(shuffle_auc)
df = data.frame(x=shuffle_auc)



fig = ggplot(data=df,aes(x=x)) + geom_histogram(alpha=0.3,colour='black',binwidth=0.1)+ geom_vline(xintercept=real_auc)


dir.create(file.path(output_folder,create_folder))
setwd(file.path(output_folder,create_folder))
pdf("shuffle_compare.pdf",width=8)
print(fig)
dev.off()



