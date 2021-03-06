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
test_type = "wilcox" # "ttest"/"wilcox"/"regression"
sig_gene = 50
p_thresh=0.001
q_thresh=0.05
p_step=0.001
q_step=0.05
p_up = 0.05
q_up = 0.2

#basic parameters for elastic net#
BS = 100
alphas = seq(0.1,1,by=0.1)

#minimum number of selected features#
feature_min = 10
freq_step = 0.05
freq = 0.8


#data preprocess#
filter_low_exp = FALSE
exp_normalize = FALSE
add_clinic = FALSE


#single cancer#



#all genes or pre-defined gene set#






if( input_type == "molecular_only" )
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
}

if( input_type == "clinical_only" )
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
  test_type = "linear_regression" # "ttest"/"wilcox"/"regression"/"linear_regression"
  sig_gene = 50
  p_thresh=0.05
  q_thresh=0.05
  p_step=0.01
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
}


