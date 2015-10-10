###parameters###
#libraries#
library(glmnet)
library(doParallel)
library(foreach)
library(pracma)
no_cores = detectCores()


#from the *th column in the info file is train/test partition
info_col = 3

##function parameters initialization##

#identify differential genes#
find_diff_genes = TRUE
test_type = "regression" # "ttest"/"wilcox"/"regression"
if(calc_cancer=="sin_cancer"){
  multi_cancer = F
}else{
  multi_cancer = T
}
sig_gene = 1
p_thresh=0.05
q_thresh=0.05
p_step=0.01
q_step=0.05
p_up = 0.1
q_up = 0.2

#basic parameters for elastic net#
BS = 100
alphas = seq(0.1,1,by=0.1)

#minimum number of selected features#
feature_min = 2
select_gene = "by_freq"
freq_step = 0.05
freq = 0.8
select_gene = "freqNweight" #should be chosen, by here for ease of use
sd = 2
sd_step = 0.1


#key parameters to output
if( select_gene == "freqNweight")
{
  key_param = data.frame(param=c("p.thresh","sig.gene","sd","freq.gene"),
                         value=rep(NA,4) )
}else{
  key_param = data.frame(param=c("p.thresh","sig.gene","freq.thresh","freq.gene"),
                         value=rep(NA,4))
}
#data preprocess#
filter_low_exp = FALSE
exp_normalize = FALSE
add_clinic = FALSE
combine_clinic = FALSE #identify markers only, predict combined with cancer type



if( input_type == "molecular_only" )
{
  
  find_diff_genes = TRUE
  test_type = "wilcox" # "ttest"/"wilcox"/"regression"
  sig_gene = 50
  p_thresh=0.05
  q_thresh=0.05
  p_step=0.01
  q_step=0.05
  p_up = 0.1
  q_up = 0.2
  
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = FALSE
  combine_clinic = FALSE
}

if( input_type == "clinical_only" )
{
  find_diff_genes = FALSE
  test_type = NULL # "ttest"/"wilcox"/"regression"
  
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = TRUE
}

if(input_type == "clinical_pred")
{
  find_diff_genes = FALSE
  test_type = NULL # "ttest"/"wilcox"/"regression"
  
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = TRUE
}

if( input_type == "clinical_molecular" )
{
  find_diff_genes = TRUE
  test_type = "wilcox" # "ttest"/"wilcox"/"regression"
  sig_gene = 50
  p_thresh=0.001
  q_thresh=0.05
  p_step=0.002
  q_step=0.05
  p_up = 0.1
  q_up = 0.2
  
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = TRUE
}

if(input_type == "half_clinical_molecular")
{
  find_diff_genes = TRUE
  test_type = "regression" # "ttest"/"wilcox"/"regression"
  sig_gene = 2
  p_thresh=0.05
  q_thresh=0.05
  p_step=0.01
  q_step=0.05
  p_up = 0.1
  q_up = 0.2
  
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = FALSE
}

if(input_type=="markers_combine_type")
{
  find_diff_genes = TRUE
  test_type = "wilcox" # "ttest"/"wilcox"/"regression"
  sig_gene = 50
  p_thresh=0.001
  q_thresh=0.05
  p_step=0.001
  q_step=0.05
  p_up = 0.05
  q_up = 0.2
  
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = FALSE
  combine_clinic = TRUE
}





setwd("functions/")
source("bootstrap_sample.R")
source("partition_data.R")
source("comb.R")
source("filter_mRNA.R")
source("exp_norm.R")
source("test_gene.R")
source("impute_NA.R")
source("cancer_dummy.R")
source("dummy_to_test.R")
source("gene_selection.R")
source("ensemble_roc.R")
source("map_rna_gene.R")
source("gene_selection_sd.R")


