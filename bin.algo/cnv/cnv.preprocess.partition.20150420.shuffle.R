##shuffle the response##
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/cnv/shuffle/")
cisplatin.info = read.table("../cisplatin.gistic2.5_fold_cv.mat.txt",header=T,sep="\t",quote="")
#View(cisplatin.info)

cancer_type = as.character(unique(cisplatin.info$cancer))

for(j in 1:100)
{
  cisplatin.alter = cisplatin.info
  for( i in 1:length(cancer_type) )
  {
    ix = which(cisplatin.info$cancer==cancer_type[i])
    curr.info = cisplatin.info[ix,]
    
    set.seed(i)
    resp = sample( as.character(curr.info$response),replace=F,size=nrow(curr.info)  )
    
    curr.info$response = resp
    
    cisplatin.alter[ix,] = curr.info
  }
  tmp_str = paste("cisplatin.gistic2.shuffle_",j,".txt",sep="")
  write.table(cisplatin.alter,tmp_str,row.names=F,col.names=T,quote=F,sep="\t")
}

