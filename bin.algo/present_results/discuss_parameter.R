#libraries#
library(ggplot2)
#functions#
plot_param <- function(df){
  p1 = ggplot(data=df,aes(x=as.factor(p_set),y=auc)) + geom_boxplot() + labs(title="AUC distribution",y="AUC",x="set p value")
  p2 = ggplot(data=df,aes(x=as.factor(p_set),y=p.thresh)) + geom_boxplot() + labs(title="Practical p value distribution",x="set p value")
  p3 = ggplot(data=df,aes(x=as.factor(p_set),y=sig.gene)) + geom_boxplot() + labs(title="Practical p value distribution",x="set p value")
  p4 = ggplot(data=df,aes(x=as.factor(p_set),y=thresh)) + geom_boxplot() + labs(title="Practical standard deviation distribution",x="set p value")
  p5 = ggplot(data=df,aes(x=as.factor(p_set),y=freq.gene)) + geom_boxplot() + labs(title="Genes/probes in genomic signature",x="set p value")
  

  #   print(p1)
  #   print(p2)
  #   print(p3)
  #   print(p4)
  #   print(p5)
  multiplot(p1, p2, p3, p4, p5, cols=2)
}

source("source_all.R")
#load data#
#desktop
# setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/cisplatin/mRNA/elastic_net/logistic.cancer/present_results/performance/")
# file_0.05 = "whole_genome.sd.0.05/all_res.99.txt"
# file_0.01 = "whole_genome.sd.0.01/all_res.100.txt"
# file_0.001 = "whole_genome.sd.0.001/all_res.99.txt"
#cluster
# args <- commandArgs(trailingOnly=TRUE)
# file_0.05 = args[1]
# file_0.01 = args[2]
# file_0.001 = args[3]
#common

##cisplatin##
#mRNA#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/cisplatin/mRNA/elastic_net/logistic.cancer/present_results/performance/")
#pan-cancer
file_0.05 = "whole_genome.sd.0.05/all_res.99.txt"
file_0.01 = "whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "whole_genome.sd.0.001/all_res.99.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#BLCA
file_0.05 = "blca.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "blca.whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "blca.whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#CESC
file_0.05 = "cesc.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "cesc.whole_genome.sd.0.01/all_res.99.txt"
file_0.001 = "cesc.whole_genome.sd.0.001/all_res.99.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#LUAD
file_0.05 = "luad.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "luad.whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "luad.whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#cnv#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/cisplatin/cnv/elanet/logistic.cancer/present_results/performance/")
#pan-cancer
file_0.05 = "whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#BLCA
#no need
file_0.05 = "blca.whole_genome.sd.0.05/all_res.71.txt"
# file_0.01 = "blca"
# file_0.001 = "blca.whole_genome.sd.0.001/"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
# dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
# dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = rep(0.05,nrow(dat_0.05))
df = cbind(dat_0.05,p_set)
plot_param(df)
#CESC
file_0.05 = "cesc.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "cesc.whole_genome.sd.0.01/all_res.99.txt"
file_0.001 = "cesc.whole_genome.sd.0.001/all_res.97.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#LUAD
file_0.05 = "luad.whole_genome.sd.0.05/all_res.99.txt"
file_0.01 = "luad.whole_genome.sd.0.01/all_res.97.txt"
file_0.001 = "luad.whole_genome.sd.0.001/all_res.94.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#methylation#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/cisplatin/methylation/elastic_net/logistic.cancer/present_results/performance/")
#pan-cancer
file_0.05 = "whole_genome.sd.0.05/all_res.97.txt"
file_0.01 = "whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "whole_genome.sd.0.001/all_res.97.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#blca
file_0.05 = "blca.whole_genome.sd.0.05/all_res.97.txt"
file_0.01 = "blca.whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "blca.whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#luad
file_0.05 = "luad.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "luad.whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "luad.whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#miRNA#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/cisplatin/miRNA/elastic_net/logistic.cancer/present_results/performance/")
#pan-cancer
file_0.05 = "whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#blca
file_0.05 = "blca.whole_genome.sd.0.05/all_res.99.txt"
file_0.01 = "blca.whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "blca.whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#cesc
file_0.05 = "cesc.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "cesc.whole_genome.sd.0.01/all_res.99.txt"
file_0.001 = "cesc.whole_genome.sd.0.001/all_res.99.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#miRNA
file_0.05 = "luad.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "luad.whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "luad.whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)



##carboplatin##
#mRNA#
#pan-cancer

#ucec

#miRNA


#methylation


#cnv#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/carboplatin/cnv/elastic_net/logistic.cancer/present_results/performance/")
#pan-cancer
file_0.05 = "whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "whole_genome.sd.0.01/all_res.99.txt"
file_0.001 = "whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)
#ucec
file_0.05 = "ucec.whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "ucec.whole_genome.sd.0.01/all_res.100.txt"
file_0.001 = "ucec.whole_genome.sd.0.001/all_res.100.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)



##paclitaxel##
#cnv#
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/paclitaxel/cnv/elastic_net/logistic.cancer/present_results/performance/")
#pan-cancer
file_0.05 = "whole_genome.sd.0.05/all_res.100.txt"
file_0.01 = "whole_genome.sd.0.01/all_res.99.txt"
file_0.001 = "whole_genome.sd.0.001/all_res.99.txt"
dat_0.05 = read.table(file_0.05,header=T,sep="\t",quote="")
dat_0.01 = read.table(file_0.01,header=T,sep="\t",quote="")
dat_0.001 = read.table(file_0.001,header=T,sep="\t",quote="")
p_set = c(rep(0.05,nrow(dat_0.05)),rep(0.01,nrow(dat_0.01)),rep(0.001,nrow(dat_0.001)))
df = cbind(rbind(dat_0.05,dat_0.01,dat_0.001),p_set)
plot_param(df)

#mRNA
setwd("C:/Users/zding/workspace/projects/drug_sensitivity/results/paclitaxel/")

#miRNA


#methylation
