#functions#
read_table <- function(x){
  return( read.table(x,sep="\t",header=T,quote="") )
}

###load data###
args <- commandArgs(trailingOnly=TRUE)
#input
input_folder = args[1] #no "/" at last

#output
output_folder = args[2]
create_folder = args[3]



dirs = list.dirs(path=input_folder,full.names=T,recursive=F)
output_names = lapply(dirs,
                      function(x){ tmp = strsplit(x,split="\\/"); 
                                   tmp=unlist(tmp); 
                                   return(tmp[length(tmp)]) } 
                      )
dir.create(file.path(output_folder,create_folder),showWarnings=FALSE)
setwd( file.path(output_folder,create_folder) )
if(length(dirs)>0){
  for( j in 1:length(dirs))
  {
    curr_folder = dirs[j]
    output_file_name = output_names[[j]]
    
    pan.files = list.files(path=curr_folder,full.names=T,pattern=param_file_pattern)
    pan.files = sort(pan.files)
    param = lapply(pan.files,read_table)

    ##get paramters##
    p.thresh = unlist( lapply(param,function(x)return(x[1,2])) )
    sig.gene = unlist( lapply(param,function(x)return(x[2,2])) )
    thresh = unlist( lapply(param,function(x)return(x[3,2])) )
    freq.gene = unlist( lapply(param,function(x)return(x[4,2])) )
    
    
    df = data.frame(p.thresh=p.thresh,sig.gene=sig.gene,
                    thresh=thresh,freq.gene = freq.gene)
    
    write.table(df,output_file_name,row.names=F,col.names=T,sep="\t",quote=F)
  }
}else{
  pan.files = list.files(path=input_folder,full.names=T,pattern=param_file_pattern)
  pan.files = sort(pan.files)
  param = lapply(pan.files,read_table)
  
  ##get paramters##
  p.thresh = unlist( lapply(param,function(x)return(x[1,2])) )
  sig.gene = unlist( lapply(param,function(x)return(x[2,2])) )
  thresh = unlist( lapply(param,function(x)return(x[3,2])) )
  freq.gene = unlist( lapply(param,function(x)return(x[4,2])) )
  
  
  df = data.frame(p.thresh=p.thresh,sig.gene=sig.gene,
                  thresh=thresh,freq.gene = freq.gene)
  
  write.table(df,"param.txt",row.names=F,col.names=T,sep="\t",quote=F)
}



