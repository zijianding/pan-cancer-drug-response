###libraries###
library(reshape2)
library(ggplot2)
library(plyr)
library(gplots)
library(NMF)


###functions###
impute_NA <- function(train_dat,gene_num)
{
  #input data: row as genes and col as samples
  #output data:row as genes and col as samples
  #usage: before glmnet function
  curr_train = train_dat[1:gene_num,]
  
  
  mean_nums = rowMeans( curr_train,na.rm=T )
  for(i in 1:nrow(train_dat))
  {
    train_na = which( is.na(train_dat[i,]) )
    
    
    if( length(train_na)>0)
    {
      curr_train[i,train_na] = mean_nums[i]
    }
    
  }
  
  train_dat[1:gene_num,] = curr_train
  return(train_dat)
  
}


###output data###
pdf("cisplatin.mRNA.litmining_marker.combine.pdf",width=14,height=7)

###load data###
##desktop##
#mRNA#
dat_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/mRNAseq/cisplatin.mRNAseq.gdac_20141206.preprocess.txt"
info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/mRNAseq/cisplatin.mRNAseq_fold_cv.mat.txt"

gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/marker_gene/by_data/marker_molecular_cancer_type.txt"
#gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/gene_sets/cgp.2012.txt"
#gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/gene_sets/gdsc.txt"

#miRNA molecular only#
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/miRNAseq/marker_gene/marker_molecular_only.txt"
# dat_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/miRNAseq/cisplatin.miRNAseq.gdac_20141206.preprocess.txt"
# info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/miRNAseq/cisplatin.miRNAseq_fold_cv.mat.txt"
##cluster##

##both##
genes_pre = read.table(gene_file,header=F,row.names=NULL,sep="\t",quote="")
cisplatin.dat = read.table(dat_file ,header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(info_file,sep="\t",header=T,quote="")
freq_thresh = 0.8

###preprocess data###
##find data##
pats = strsplit(colnames(cisplatin.dat),"\\.")
pats_cancer = unlist(lapply(pats,function(x){return(x[[1]][1])}))

# dat_genes = strsplit(rownames(cisplatin.dat),"\\|")
# dat_genes = unlist(lapply(dat_genes,function(x){return(x[[1]][1])}))
# genes = genes_pre[,1][genes_pre[,1] %in% dat_genes]

genes = genes_pre[,1][genes_pre[,2]>=freq_thresh]
cancer_ix = which(genes %in% unique(pats_cancer))
if(length(cancer_ix)>0)
{
  cancer_mat = matrix(NA,nrow=length(cancer_ix),ncol=length(pats_cancer))
  for(i in 1:length(cancer_ix))
  {
    cancer_mat[i,which(genes[cancer_ix[i]]==pats_cancer)] = 1
    cancer_mat[i,which(genes[cancer_ix[i]]!=pats_cancer)] = 0
  }
  rownames(cancer_mat) = genes[cancer_ix]
  colnames(cancer_mat) = colnames(cisplatin.dat)
}
gene_ix = match(genes,rownames(cisplatin.dat))
#gene_ix = match(genes,dat_genes)
gene_ix = gene_ix[!is.na(gene_ix)]


pre_dat = cisplatin.dat[gene_ix,]
if(length(cancer_ix)>0)
{
  pre_dat = rbind(pre_dat,cancer_mat)
}


resp = as.character(cisplatin.info$response[match(colnames(cisplatin.dat),cisplatin.info$patient)])
resp[resp=="insensitive"] = "non-responder"
resp[resp=="sensitive"] = "responder"

resp_sort = sort(resp,index=T)
resp = resp_sort$x
pre_dat = pre_dat[,resp_sort$ix]
pats_cancer = pats_cancer[resp_sort$ix]


final_dat = impute_NA(pre_dat,nrow(pre_dat))


###draw heatmap###
#heatmap#
colAnno = data.frame(response=resp,cancer=pats_cancer)
row_lab.list = strsplit(rownames(final_dat),"\\|")
row_lab = unlist(lapply(row_lab.list,function(x){return(x[[1]][1])}))
aheatmap(final_dat,annCol=colAnno,hclustfun='average', Colv = NA, labRow=row_lab,
         distfun=function(x) as.dist((1-cor(t(x)))/2),scale="row",legend=F)
#in each cancer heatmap
cancer_type = unique(pats_cancer)
i = 4
ix = which(pats_cancer==cancer_type[i])
colAnno = data.frame(response=resp[ix],cancer=pats_cancer[ix])
aheatmap(final_dat[,ix],annCol=colAnno,hclustfun='average', Colv = NA, labRow=row_lab,
         distfun=function(x) as.dist((1-cor(t(x)))/2),scale="row",legend=F)


##analyse some genes##
df = data.frame(pats_cancer,resp,
                t(final_dat[match("PRRG1|5638",rownames(final_dat)),]))
colnames(df) = c("cancer","response","NQO1")
ggplot(df,aes(y=NQO1,x=cancer)) + geom_boxplot(aes(fill=response))
#geom_violin(aes(fill=response),alpha=1,width=1) 



#+ geom_jitter(aes(fill=response))+ +




#cc = c(rep("#FF0000FF",length(sen_ix)),rep("#CCFF00FF",times=ncol(pre_dat)-length(sen_ix)))
# cc = vector(mode="character",length=length(resp))
# cc[resp=="non-responder"] = "#CCFF00FF"
# cc[resp=="responder"] = "#FF0000FF"
# exp_heatmap = heatmap.2(data.matrix(final_dat),Rowv=T,Colv=T,scale="none",key=T,,
#                         ColSideColors=cc,col=cm.colors(255),labCol = NA,
#                         hclust=function(x) hclust(x,method="average"),
#                         distfun=function(x) as.dist((1-cor(t(x)))/2))



##NMF##

#non-negative matrix factorization#
final_dat = t(2^final_dat)
res.estimate = nmf(t(final_dat),2:10,nrun=10,seed=123456)
plot(res.estimate)
res = nmf(t(final_dat),3)
basismap(res)
coefmap(res,annCol=colAnno)
consensusmap(res,annCol=colAnno)

dev.off()



