###libraries###
library(reshape2)
library(ggplot2)
library(plyr)
library(gplots)
library(NMF)



###load data###
##desktop##
#mRNA#
data_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/data/cisplatin.mRNAseq.gdsc.preprocess.txt"
info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/data/cisplatin.mRNAseq_fold_cv.gdsc.txt"
cisplatin.dat = read.table(data_file ,header=T,row.names=1,sep="\t",quote="")
cisplatin.info = read.table(info_file,sep="\t",header=T,quote="")
##add cancer type##
#from method
#regression
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/results/marker.no_preselection.txt"
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/results/marker.clinical_molecular.0.001.txt"
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/results/marker.clinical_molecular.0.01.txt"
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/results/marker.clinical_molecular.0.05.txt"
#spearman cor
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/results/marker.clinical_molecular.0.001.spearman_cor.txt"
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/results/marker.clinical_molecular.0.01.spearman_cor.txt"
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/pan-cancer-drug-response/bin.algo/coding_correct/GDSC/results/marker.clinical_molecular.0.05.spearman_cor.txt"
#from literature
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/gene_sets/cgp.2012.txt"
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/gene_sets/gdsc.txt"
##no cancer type##
#miRNA molecular only#
# gene_file = "C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/miRNAseq/marker_gene/marker_molecular_only.txt"
# dat_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/miRNAseq/cisplatin.miRNAseq.gdac_20141206.preprocess.txt"
# info_file = "C:/Users/zding/workspace/projects/drug_sensitivity/data/omics.drug_centric/miRNAseq/cisplatin.miRNAseq_fold_cv.mat.txt"
##cluster##

##both##
genes_pre = read.table(gene_file,header=F,row.names=NULL,sep="\t",quote="")
freq_thresh = 0.8

pats = as.character(colnames(cisplatin.dat))
pats_cancer = as.character(cisplatin.info$cancer[match(pats,cisplatin.info$patient)])

gene_freq = genes_pre[,2]
gene_thresh = mean(gene_freq) + 2*sd(gene_freq)


genes = as.character(genes_pre[,1][genes_pre[,2]>=freq_thresh])
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
gene_ix = gene_ix[!is.na(gene_ix)]


pre_dat = cisplatin.dat[gene_ix,]
if(length(cancer_ix)>0)
{
  pre_dat = rbind(pre_dat,cancer_mat)
}


resp = cisplatin.info$response[match(colnames(cisplatin.dat),cisplatin.info$patient)]
resp = as.vector(resp)

resp_sort = sort(resp,index=T)
resp = resp_sort$x
pre_dat = pre_dat[,resp_sort$ix]
pats_cancer = pats_cancer[resp_sort$ix]

final_dat = pre_dat


###draw heatmap###
###output data###
tmp_str = paste(gene_file,"pdf",sep=".")
pdf(tmp_str,width=14,height=7)

#heatmap#
colAnno = data.frame(IC50=resp,cancer=pats_cancer)
row_lab = rownames(final_dat)
aheatmap(final_dat,annCol=colAnno,hclustfun='average', Colv = NA, labRow=row_lab,
         distfun=function(x) as.dist((1-cor(t(x)))/2),legend=F,scale="row")
aheatmap(final_dat,annCol=colAnno,hclustfun='average', labRow=row_lab,
         distfun=function(x) as.dist((1-cor(t(x)))/2),legend=F,scale="row")
#try the most sensitive and most insensitive each 20 samples
ix_1 = seq(1,20,by=1)
ix_2 = seq(length(pats_cancer)-19,length(pats_cancer),1)
ix = c(ix_1,ix_2)
aheatmap(final_dat[,ix],annCol=colAnno[ix,],hclustfun='average', Colv = NA, labRow=row_lab,
         distfun=function(x) as.dist((1-cor(t(x)))/2),legend=F,scale="row")
aheatmap(final_dat[,ix],annCol=colAnno[ix,],hclustfun='average', labRow=row_lab,
         distfun=function(x) as.dist((1-cor(t(x)))/2),legend=F,scale="row")
dev.off()
#in each cancer heatmap
# cancer_type = unique(pats_cancer)
# i = 4
# ix = which(pats_cancer==cancer_type[i])
# colAnno = data.frame(response=resp[ix],cancer=pats_cancer[ix])
# aheatmap(final_dat[,ix],annCol=colAnno,hclustfun='average', Colv = NA, labRow=row_lab,
#          distfun=function(x) as.dist((1-cor(t(x)))/2),legend=F)


##analyse some genes##
# gene = "MAGEA4"
# df = data.frame(pats_cancer,resp,
#                 t(final_dat[match(gene,rownames(final_dat)),]))
# colnames(df) = c("cancer","response",gene)
# ggplot(df,aes(y=MAGEA4,x=cancer)) + geom_boxplot(aes(fill=response))
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
# final_dat = t(2^final_dat)
# res.estimate = nmf(t(final_dat),2:10,nrun=10,seed=123456)
# plot(res.estimate)
# res = nmf(t(final_dat),3)
# basismap(res)
# coefmap(res,annCol=colAnno)
# consensusmap(res,annCol=colAnno)





