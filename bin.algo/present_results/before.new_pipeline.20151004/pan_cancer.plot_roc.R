###mRNA###
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/bin.algo/mRNASeq/elastic_net.logistic")
molecular = read.table("molecular_only.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("combine_molecular_cancer_type.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("cancer_type_only.txt",header=F,row.names=NULL,sep="\t")

tiff("molecular_cancer_comparision.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("topleft",lty=c(1,1,1,2),col=c("red","green","black","gray") ,
       legend=c("molcular only","cancer type only","both","random"))
title("Pan-cancer Classification of Drug response","mRNA-seq")
dev.off()

#full data elastic net#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/performance/by_data/")
molecular = read.table("molecular_only.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("combine_molecular_cancer_type.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("cancer_type_only.txt",header=F,row.names=NULL,sep="\t")
full_mc = read.table("full_dat.molecular_cancer.txt",header=F,row.names=NULL,sep="\t")

tiff("full_data.molecular_cancer_comparision.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(full_mc[,2],full_mc[,1],col="blue",lty=1)
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("topleft",lty=c(1,1,1,1,2),col=c("red","green","black","blue","gray") ,
       legend=c("molcular only","cancer type only","both","full data","random"))
title("Pan-cancer Classification of Drug response","mRNA-seq")
dev.off()


###miRNA###
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/omics_feature/miRNAseq/elastic_net")
molecular = read.table("molecular_only.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("combine_molecular_cancer_type.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("cancer_type_only.txt",header=F,row.names=NULL,sep="\t")

tiff("molecular_cancer_comparision.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("topleft",lty=c(1,1,1,2),col=c("red","green","black","gray") ,
       legend=c("molcular only","cancer type only","both","random"))
title("Pan-cancer Classification of Drug response","miRNA-seq")
dev.off()



###focal CNV mapped to genes###
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/bin.algo/cnv/elastic_net.logistic/")
molecular = read.table("pan_roc.cnv_only.elanet.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("pan_roc.cnv.elanet.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("pan_roc.cancer_only.elanet.txt",header=F,row.names=NULL,sep="\t")
tiff("molecular_cancer_comparision.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("topleft",lty=c(1,1,1,2),col=c("red","green","black","gray") ,
       legend=c("molcular only","cancer type only","both","random"))
title("Pan-cancer Classification of Drug response","Focal CNV mapped to genes")
dev.off()


###methylation###
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/methylation/")
molecular = read.table("molecular_only.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("combine_molecular_cancer_type.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("cancer_type_only.txt",header=F,row.names=NULL,sep="\t")

tiff("molecular_cancer_comparision.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("topleft",lty=c(1,1,1,2),col=c("red","green","black","gray") ,
       legend=c("molcular only","cancer type only","both","random"))
title("Pan-cancer Classification of Drug response","DNA methylation")
dev.off()


###mRNA known gene sets###
##han re zheng##
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/known_gene_sets/hanrezheng/")
molecular = read.table("hanrezheng_molecular_only.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("hanrezheng_both.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/bin.algo/mRNASeq/elastic_net.logistic/cancer_type_only.txt",
                    header=F,row.names=NULL,sep="\t")
data_sole = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/bin.algo/mRNASeq/elastic_net.logistic/molecular_only.txt",
                       header=F,row.names=NULL,sep="\t")
tiff("molecular_cancer_data_sole.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(data_sole[,2],data_sole[,1],col="blue",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("bottomright",lty=c(1,1,1,1,2),col=c("red","green","black","blue","gray") ,
       legend=c("HanReZheng Genes","cancer type only","HanReZheng and cancer","molecular only","random"))
title("Pan-cancer Classification of Drug response","mRNA-seq")
dev.off()


##cytoscape cisplatin genes##
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/known_gene_sets/cytoscape_cisplatin_genes/")
molecular = read.table("cytoscape_cisplatin_gene_molecular_only.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("cytoscape_cisplatin_gene_both.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/bin.algo/mRNASeq/elastic_net.logistic/cancer_type_only.txt",
                    header=F,row.names=NULL,sep="\t")
data_sole = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/bin.algo/mRNASeq/elastic_net.logistic/molecular_only.txt",
                       header=F,row.names=NULL,sep="\t")
tiff("molecular_cancer_data_sole.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(data_sole[,2],data_sole[,1],col="blue",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("bottomright",lty=c(1,1,1,1,2),col=c("red","green","black","blue","gray") ,
       legend=c("Literature Mining Cisplatin related Genes","cancer type only","Mining Genes and cancer","molecular only","random"))
title("Pan-cancer Classification of Drug response","mRNA-seq")
dev.off()


##cpg 2012 cisplatin mRNA markers##
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/performance/cgp2012")
molecular = read.table("cgp2012_molecular.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("cgp2012_both.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/performance/by_data/cancer_type_only.txt",
                    header=F,row.names=NULL,sep="\t")
data_sole = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/performance/by_data/molecular_only.txt",
                       header=F,row.names=NULL,sep="\t")
tiff("molecular_cancer_data_sole.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(data_sole[,2],data_sole[,1],col="blue",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("bottomright",lty=c(1,1,1,1,2),col=c("red","green","black","blue","gray") ,
       legend=c("CGP2012 genes","cancer type only","CGP2012 Genes and cancer","molecular only","random"))
title("Pan-cancer Classification of Drug response","mRNA-seq")
dev.off()


###gdsc cisplatin mRNA markers###
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/performance/gdsc/")
molecular = read.table("gdsc_molelcular.txt",header=F,row.names=NULL,sep="\t")
molecular_cancer = read.table("gdsc_both.txt",header=F,row.names=NULL,sep="\t")
cancer = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/performance/by_data/cancer_type_only.txt",
                    header=F,row.names=NULL,sep="\t")
data_sole = read.table("C:/Users/zding/workspace/projects/drug_sensitivity/results/present_results/mRNA/performance/by_data/molecular_only.txt",
                       header=F,row.names=NULL,sep="\t")
tiff("molecular_cancer_data_sole.tiff")
plot(molecular[,2],molecular[,1],col="red",lty=1,"l",
     xlab="False Positive Rate",ylab="True Positive Rate")
lines(cancer[,2],cancer[,1],col="green",lty=1,"l")
lines(molecular_cancer[,2],molecular_cancer[,1],col="black",lty=1,"l")
lines(data_sole[,2],data_sole[,1],col="blue",lty=1,"l")
lines(seq(0,1,by=0.1),seq(0,1,by=0.1),lty=2,"l",col="gray")
legend("bottomright",lty=c(1,1,1,1,2),col=c("red","green","black","blue","gray") ,
       legend=c("GDSC genes","cancer type only","GDSC Genes and cancer","molecular only","random"))
title("Pan-cancer Classification of Drug response","mRNA-seq")
dev.off()
