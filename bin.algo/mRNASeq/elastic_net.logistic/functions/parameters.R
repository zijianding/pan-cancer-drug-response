###parameters###
#basic parameters for elastic net#
BS = 100
alphas = seq(0.1,1,by=0.1)

#minimum number of selected features#
feature_min = 10
freq_step = 0.05
freq = 0.8

#from the *th column in the info file is train/test partition
info_col = 3

#initialization
find_diff_genes = TRUE
diff_type = "ttest" # "ttest"/"wilcox"/"regression"
filter_low_exp = FALSE
exp_normalize = FALSE
add_clinic = FALSE



if( input_type == "molecular_only" )
{
  find_diff_genes = FALSE
  diff_type = "wilcox" # "ttest"/"wilcox"/"regression"
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = FALSE
}

if( input_type == "clinical_only" )
{
  find_diff_genes = FALSE
  diff_type = NULL # "ttest"/"wilcox"/"regression"
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = TRUE
}

if( input_type == "clinical_molecular" )
{
  find_diff_genes = TRUE
  diff_type = "wilcox" # "ttest"/"wilcox"/"regression"
  filter_low_exp = FALSE
  exp_normalize = FALSE
  add_clinic = TRUE
}


