ef_method_pipeline <- function(counts_list, reference_standard, subject_df, bonfer_correction = 1){
#Description: It calculates the natural log odds-ratio using Gart's method and their associated variances. It then runs the 
              #EF-method and produces results for all GO terms. It also produces p-values which have their FDR adjustment based
              #only on the subset of 4 listed in the reference standard. 
#Input: counts_list (list of matrices) - the counts for the contingency table for the GO terms. The 1st column must be genes with
                                        #altered expression inside the GO term, 2nd must be genes with unaltered expression in
                                        #GO term, 3rd must be genes with altered expression outside GO term, and 4th must be 
                                        #genes with unaltered expression outside the GO term. 
       #reference_standard (data.frame) - contains the GO terms in the reference standard and whether or not they are 
                                          #'dysregulated' or in this case just have altered expression. 
       #subject_df (data.frame) - describes each subject. Must have a variable called 'Subject' which provides subject unique 
                                  #name and variable 'Mutation' which provides which group the subject belongs to. 
       #bonfer_correction (double) - The Bonferroni correction for the nominal p-values that you want to apply. It must be 
                                    #positive and greater than or equal to 1. 
#Output: output_list (list of dataframe/lists) - contains all the info produced from the pipeline including diagnostics. 1st object
                                                #is called 'subset_res' and contains the info related to the subset of GO terms
                                                #in the reference standard. Here the rows are named for GO terms, and the last
                                                #column is the p-values adjusted via FDR using only the 4 present in the 
                                                #reference standard. 2nd object is called 'ef_res' which provides the summary of
                                                #results for all the GO terms. 3rd object 'ef_info' contains all the diagnostic
                                                #info for the EF-method run. 

# Calculate the log-OR and the variance of the log-ORs
  # Function: calc_lnOR_gart
  subject_lnors <- lapply(counts_list, calc_lnOR_gart)
  
# Combine the log-OR info across subjects and then do the same for their associated variances
lnor_mat <- sapply(subject_lnors, function(x){ return(x[, 1]) })
var_lnor_mat <- sapply(subject_lnors, function(x){ return(x[, 2]) })

# If the data only has one GO terms use a different version of the function designed for vectors
if(!is.matrix(lnor_mat)){
  ef_info <- run_ef_method_single_gs(lnor_mat = lnor_mat, var_lnor_mat = var_lnor_mat, subject_df = subject_df, 
                                      bonfer_correction = bonfer_correction)
  ef_res <- ef_info$summary
  rownames(ef_res) <- rownames(subject_lnors[[1]])
} else{
  # Run the EF-method
  ef_info <- run_ef_method(lnor_mat = lnor_mat, var_lnor_mat = var_lnor_mat, subject_df = subject_df, 
                           bonfer_correction = bonfer_correction)
  ef_res <- ef_info$summary
}

# Subset the GO terms to those in the reference standard
subset_res <- ef_res[rownames(ef_res) %in% reference_standard$goid, ]
fdr_vals_based_subset <- p.adjust(subset_res$bonfer_p_value, method = 'BH')
fdr_vals_based_subset[!subset_res$enriched] <- 1 # P-values adjusted based on only the subset of 4
subset_res <- cbind(subset_res, fdr_vals_based_subset)

# Combine all the results and output
output_list <- list(subset_res = subset_res, ef_res = ef_res, ef_info = ef_info)
return(output_list)
}

calc_lnOR_gart <- function(count_matrix){
#Description: The function takes a count matrix and calculates the natural log odds-ratios and their associated variances of the 
              #contingency table:         Altered     Unaltered
                            # Inside GO:   c11          c12
                            #Outside GO:   c21          c22
#Input: count_matrix (matrix of integers) - the counts for the contingency table for each GO term. The rownames must be GO terms and
                                            #the 1st column must be the count of genes with altered expression inside the GO term,
                                            #the 2nd column must be the count of genes with unaltered expression inside the GO term,
                                            #the 3rd column must be the count of genes with altered expression outside the GO term,
                                            #and the 4th column must be the count of genes with unaltered expression outside the GO
                                            #term. 
#Output: lnor_info (matrix) - the log odds-ratios and their associated variances, which were calculated using Gart's method. The
                              #rows are named by GO terms and correspond to such. The 1st column is named "lnOR" and contains the
                              #calculated log odds-ratios. The 2nd column is named "var_lnOR" and contains the variance of the 
                              #associated log odds-ratios. 
  
# Calculate the lnOR, with adding 0.5 to each count
lnORs <- log( ((count_matrix[, 1] + .5) * (count_matrix[, 4] + .5)) / ((count_matrix[, 2] + .5) * (count_matrix[, 3] + .5)) )

# Calculate the variance of the lnORs
var_lnORs <- (1 / (count_matrix[, 1] + .5)) + (1 / (count_matrix[, 2] + .5)) + (1 / (count_matrix[, 3] + .5)) + 
              (1 / (count_matrix[, 4] + .5))

# Combine and name the output
lnor_info <- cbind(lnORs, var_lnORs)
colnames(lnor_info) <- c('lnOR', 'var_lnOR')
row.names(lnor_info) <- rownames(count_matrix)
return(lnor_info)
}

# Run short EF-method
run_ef_method <- function(lnor_mat, var_lnor_mat, subject_df, bonfer_correction = 1){
#Description: It runs the EF-method and calculates the p-values and fdr-adjusted p-values. It then adjusts the fdr values based
              #on whether or not the GO term is enriched. It can also Bonferroni correct the p-values just in case the analyses
              #are run for multiple directions.
#Input: lnor_mat (matrix of doubles) - contains all the natural log odds-ratios of the subjects for each of the GO terms. The 
                                      #rows should be named by GO terms and the columns should be named by subject. 
       #var_lnor_mat (matrix of doubles) - contains all the variances of the natural log odds-ratios of the subjects for each
                                          #of the GO terms. The rows should be named by GO terms and the columns should be named
                                          #by subject. 
       #subject_df (dataframe of strings) - lists the info for each of the subjects. Must contain the variable "Subject" which 
                                            #provides the subject specific names and the variable "Mutation" which provides the
                                            #group to which the subject belongs. 'Mutation' variable must be a factor and have
                                            #levels otherwise it won't work.
       #bonfer_correction (double) - The Bonferroni correction for the nominal p-values that you want to apply. It must be 
                                    #positive and greater than or equal to 1. 
#Output: ef_output (list of dataframes and list) - 1st object is called "summary" and shows the results for the EF-method. The
                                                  #columns of the data.frame are "statistic" which provides the EF-method statistic
                                                  #value for the GO term, "p_value" which provides the nominal p-values, 'fdr_value'
                                                  #which provides the FDR adjusted p-values, 'adjusted_fdr_values' which provides 
                                                  #the FDR adjusted p-values with non-enriched GO terms set to 1, and 'enriched'
                                                  #which indicates whether or not the GO term was enriched or not. 2nd object is 
                                                  #called 'stat_diag' and it contains all the group summary info. It has column
                                                  #'mean_lnor1' which provides the mean lnOR for the 1st group, 'mean_lnor2' which
                                                  #provides the mean lnOR for the 2nd group, 'var_mean_lnor1' which provides the
                                                  #variance of the mean lnOR for group 1, and 'var_mean_lnor2' which provides 
                                                  #the variance of the mean lnOR for group 2. 3rd object is called 'group_vals'
                                                  #and it provides the lnORs and variance of the subjects for within each group. 
                                                  #The 'lnors1' and 'lnors2' provide the calculated lnOR's for groups 1 and 2 and
                                                  #'var_lnors1' and 'var_lnors2' provides the the variances of the lnORs for 
                                                  #groups 1 and 2. 
  
# Find the subjects which belong to each group
poss_levels <- levels(subject_df$Mutation)
group1 <- subject_df$Subject[ subject_df$Mutation == poss_levels[1]] 
group2 <- subject_df$Subject[ subject_df$Mutation == poss_levels[2]]

# Put together the lnORs and the variances associated with each group 
  # For group 1
  lnors1 <- lnor_mat[, colnames(lnor_mat) %in% group1]
  var_lnors1 <- var_lnor_mat[, colnames(var_lnor_mat) %in% group1]
  
  # For group 2
  lnors2 <- lnor_mat[, colnames(lnor_mat) %in% group2]
  var_lnors2 <- var_lnor_mat[, colnames(var_lnor_mat) %in% group2]
  
# Calculate the average lnOR for each group
mean_lnor1 <- rowMeans(lnors1)
mean_lnor2 <- rowMeans(lnors2)

# Calculate if the GO term is enriched
is_enriched <- !((mean_lnor1 <= 0) & (mean_lnor2 <= 0))
  
# Calculate the variance of the lnOR for each group
var_mean_lnor1 <- rowSums(var_lnors1) / (length(group1))^2
var_mean_lnor2 <- rowSums(var_lnors2) / (length(group2))^2

# Calculate the EF-method statistic
  # Calculate the numerator
  numerator <- mean_lnor1 - mean_lnor2
  
  # Calculate the denominator
  denominator <- sqrt(var_mean_lnor1 + var_mean_lnor2)
  
  # Put together the statistic
  ef_stat <- numerator / denominator 
  
  # Find out the p-value and adjust it for FDR
  p_value <- 2 * pnorm(abs(ef_stat), lower.tail = FALSE)
  bonfer_p_value <- sapply(p_value, function(x){ 
    return(min(x * bonfer_correction, 1)) 
  })
  fdr_values <- p.adjust(bonfer_p_value, method = 'BH')
  
  # Set the p-values to 1 if the GO term was not enriched
  adjusted_fdr_values <- fdr_values
  adjusted_fdr_values[!is_enriched] <- 1
  
# Combine all the info and output
output_summary <- data.frame(statistic = ef_stat, p_value = p_value, bonfer_p_value = bonfer_p_value, fdr_value = fdr_values, 
                             adjusted_fdr_values = adjusted_fdr_values, enriched = is_enriched)
stat_diagnostics <- data.frame(mean_lnor1 = mean_lnor1, mean_lnor2 = mean_lnor2, var_mean_lnor1, var_mean_lnor2)
group_values <- list(lnors1 = lnors1, lnors2 = lnors2, var_lnors1 = var_lnors1, var_lnors2 = var_lnors2)
ef_output <- list(summary = output_summary, stat_diag = stat_diagnostics, group_vals = group_values)
return(ef_output)
}

run_ef_method_single_gs <- function(lnor_mat, var_lnor_mat, subject_df, bonfer_correction = 1){
#Description: It runs Inter-N-of-1 for the case where there is only one gene set which you wish to test. 
#Input: lnor_mat (vector of doubles) - contains all the natural log odds-ratios of the subjects for the single GO term. The 
                                      #objects should be named by subject. 
       #var_lnor_mat (vector of doubles) - contains all the variances of the natural log odds-ratios of the subjects for the 
                                          #single GO term. The objects should be named by subject.

       #subject_df (dataframe of strings) - lists the info for each of the subjects. Must contain the variable "Subject" which 
                                           #provides the subject specific names and the variable "Mutation" which provides the
                                           #group to which the subject belongs. 'Mutation' variable must be a factor and have
                                           #levels otherwise it won't work.
       #bonfer_correction (double) - The Bonferroni correction for the nominal p-values that you want to apply. It must be 
                                     #positive and greater than or equal to 1. 
#Output: ef_output (list of dataframes and list) - 1st object is called "summary" and shows the results for the EF-method. The
                                                  #columns of the data.frame are "statistic" which provides the EF-method statistic
                                                  #value for the GO term, "p_value" which provides the nominal p-values, 'fdr_value'
                                                  #which provides the FDR adjusted p-values, 'adjusted_fdr_values' which provides 
                                                  #the FDR adjusted p-values with non-enriched GO terms set to 1, and 'enriched'
                                                  #which indicates whether or not the GO term was enriched or not. 2nd object is 
                                                  #called 'stat_diag' and it contains all the group summary info. It has column
                                                  #'mean_lnor1' which provides the mean lnOR for the 1st group, 'mean_lnor2' which
                                                  #provides the mean lnOR for the 2nd group, 'var_mean_lnor1' which provides the
                                                  #variance of the mean lnOR for group 1, and 'var_mean_lnor2' which provides 
                                                  #the variance of the mean lnOR for group 2. 3rd object is called 'group_vals'
                                                  #and it provides the lnORs and variance of the subjects for within each group. 
                                                  #The 'lnors1' and 'lnors2' provide the calculated lnOR's for groups 1 and 2 and
                                                  #'var_lnors1' and 'var_lnors2' provides the the variances of the lnORs for 
                                                  #groups 1 and 2. 
  
# Find the subjects which belong to each group
poss_levels <- levels(subject_df$Mutation)
group1 <- subject_df$Subject[ subject_df$Mutation == poss_levels[1]] 
group2 <- subject_df$Subject[ subject_df$Mutation == poss_levels[2]]
  
# Put together the lnORs and the variances associated with each group 
  # For group 1
  lnors1 <- lnor_mat[names(lnor_mat) %in% group1]
  var_lnors1 <- var_lnor_mat[names(var_lnor_mat) %in% group1]
  
  # For group 2
  lnors2 <- lnor_mat[names(lnor_mat) %in% group2]
  var_lnors2 <- var_lnor_mat[names(var_lnor_mat) %in% group2]
  
  # Calculate the average lnOR for each group
  mean_lnor1 <- mean(lnors1)
  mean_lnor2 <- mean(lnors2)
  
  # Calculate if the GO term is enriched in at least one cohort
  is_enriched <- !((mean_lnor1 <= 0) & (mean_lnor2 <= 0))
  
  # Calculate the variance of the lnOR for each group
  var_mean_lnor1 <- sum(var_lnors1) / (length(group1))^2
  var_mean_lnor2 <- sum(var_lnors2) / (length(group2))^2
  
# Calculate the EF-method statistic
  # Calculate the numerator
  numerator <- mean_lnor1 - mean_lnor2
  
  # Calculate the denominator
  denominator <- sqrt(var_mean_lnor1 + var_mean_lnor2)
  
  # Put together the statistic
  ef_stat <- numerator / denominator 
  
# Find out the p-value and adjust it for FDR
p_value <- 2 * pnorm(abs(ef_stat), lower.tail = FALSE)
bonfer_p_value <- sapply(p_value, function(x){ 
  return(min(x * bonfer_correction, 1)) 
})
fdr_values <- p.adjust(bonfer_p_value, method = 'BH')
  
# Set the p-values to 1 if the GO term was not enriched
adjusted_fdr_values <- fdr_values
adjusted_fdr_values[!is_enriched] <- 1
  
# Combine all the info and output
output_summary <- data.frame(statistic = ef_stat, p_value = p_value, bonfer_p_value = bonfer_p_value, fdr_value = fdr_values, 
                             adjusted_fdr_values = adjusted_fdr_values, enriched = is_enriched)
stat_diagnostics <- data.frame(mean_lnor1 = mean_lnor1, mean_lnor2 = mean_lnor2, var_mean_lnor1, var_mean_lnor2)
group_values <- list(lnors1 = lnors1, lnors2 = lnors2, var_lnors1 = var_lnors1, var_lnors2 = var_lnors2)
ef_output <- list(summary = output_summary, stat_diag = stat_diagnostics, group_vals = group_values)
return(ef_output)
}

