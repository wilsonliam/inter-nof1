# Set the working directory
setwd('C:/Users/daberasturi/Dropbox/Simulation 2020/ISMB 2022 Manuscript Direction/Code Example/Website Code Example/')

# Read in the necessary libraries (these will have to be installed to perform code)
library(mixtools)
library(doParallel)
library(GO.db)
library(data.table)

# Load necessary inputs
load('./ex_inputs.rdata')

# Load example outputs to compare
load('./ex_outputs.rdata')

# Source the necessary functions found in sourced_functions
source('./MixEnrich.R')
source('./runMixEnrich.R')
source('./differential_expression_direction.R')
source('./ef_method_pipeline.R')

# Perform N-of-1 MixEnrich on the samples and obtain the contingency table counts for the different directions
contingency_table_info <- obtain_mixenrich_counts(log2fcs = log2fcs,
                                                  subj_names = colnames(log2fcs),
                                                  poster_prob_cutoff = .99,
                                                  the_alternative = 'two.sided',
                                                  anno_file = GO_anno,
                                                  abs_log2fc_cutoff = log2(1.2), num_core = 4)

# Select the appropriate contingency table counts for each direction
bidirect_subj_cont_tables <- contingency_table_info$mixenrich_counts
up_subj_cont_tables <- contingency_table_info$mixenrich_counts_up
down_subj_cont_tables <- contingency_table_info$mixenrich_counts_down

# If you are interested only in a subset of gene sets set up a reference standard containing only them
selected_go_terms <- unique(GO_anno$path_id)[1:100]
reference_standard <- GO_anno[GO_anno$path_id %in% selected_go_terms, ]
colnames(reference_standard)[1] <- c('goid')

# Perform the bidirectional Inter-N-of-1 and the up directional Inter-N-of-1 and down directional Inter-N-of-1
  # Run Bidirectional Inter-N-of-1 with no bonferroni correction of p-values
  bidirect_inter_nof1_res <- ef_method_pipeline(counts_list = bidirect_subj_cont_tables,
                                                reference_standard = reference_standard,
                                                subject_df = subj_pheno_file,
                                                bonfer_correction = 1)

    # Note: if testing bidirectional, up, and down directional all at same time then set bonfer_correction = 3
           #if testing only two directions: say up and down at same time set bonfer_correction = 2
           #if testing only 1 direction set bonfer_correction = 1

  # Run up directional Inter-N-of-1 with Bonferroni correction to account for also running Down Directional Inter-N-of-1
  up_inter_nof1_res <- ef_method_pipeline(counts_list = up_subj_cont_tables,
                               reference_standard = reference_standard,
                               subject_df = subj_pheno_file,
                               bonfer_correction = 2)

  # Run down directional Inter-N-of-1 with Bonferroni correction to account for also running Up Directional Inter-N-of-1
  down_inter_nof1_res <- ef_method_pipeline(counts_list = down_subj_cont_tables,
                                 reference_standard = reference_standard,
                                 subject_df = subj_pheno_file,
                                 bonfer_correction = 2)

# Examine the results
  # To see results for all gene sets in annotation file look in ef_res object
  head(bidirect_inter_nof1_res$ef_res)
  head(up_inter_nof1_res$ef_res)
  head(down_inter_nof1_res$ef_res)

  # To see results for subset of gene sets selected in 'reference_standard' input, look in subset_res object
  head(bidirect_inter_nof1_res$subset_res)
  head(up_inter_nof1_res$subset_res)
  head(down_inter_nof1_res$subset_res)

  # To look into inner workings of Inter-N-of-1 look into ef_info (Useful for data inspection and making plots)
  str(bidirect_inter_nof1_res$ef_info)
  str(up_inter_nof1_res$ef_info)
  str(down_inter_nof1_res$ef_info)

