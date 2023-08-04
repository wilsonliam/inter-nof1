add_side_regulation_info <- function(mix_deg_info, subj_log2fcs, genes_in_go_list){
#Description: It determines which DEGs are up and down for mixenrich, adds this info to the gene level output of mixenrich and
              #creates contingency tables based on whether the genes are upregulated or downregulated and belong in said GO terms 
#Input: mix_deg_info (dataframe) - Contains all the gene level output from MixEnrich. Must have a variable called 'symbol' which
                                  #lists the gene names and 'DEstatus' which is a binary indicator that provides whether or not
                                  #MixEnrich found the gene to be differentially expressed (1 = DEG, 0 = not DEG).
       #subj_log2fcs (vector of doubles) - provides the log2fcs of all the genes for a subject. Must have values named by their 
                                          #gene names and gene names must be same ones used for mix_deg_info. 
       #genes_in_go_list (list of vectors of strings) - provides which genes belong to which GO terms. Each object corresponds
                                                       #to a unique GO term. The names of the genes must match those used within
                                                       #mix_deg_info. 
#Output: subj_direct_info (list of dataframe/matrices) - the first object is called 'all_deg_info' and is a dataframe which 
                                                        #provides the genes, their posterior probability cutoff, their DE status,
                                                        #their upregulation DEG status, their down regulated DEG status, and
                                                        #their log2FCs. The 2nd object is called 'mix_up_counts' and provides
                                                        #the contingency table counts where the positive status is based on 
                                                        #upregulated DE Status. 3rd object is called 'mix_down_counts' and
                                                        #provides the contingency table counts where the positive status is 
                                                        #based on downregulated DE status. 2nd and 3rd object have rows based on
                                                        #GO terms and columns which provide whether the gene is positve in GO,
                                                        #negative in GO, positive outside GO, negative outside GO. Columns are 
                                                        #named in the same manner as mixenrich_counts though. 
  
# Separate out all the DEGs which were up and down
degs <- mix_deg_info[ mix_deg_info$DEstatus == 1, ]
up_genes   <- names(subj_log2fcs)[subj_log2fcs > 0]
down_genes <- names(subj_log2fcs)[subj_log2fcs <= 0]
up_degs <- degs$symbol[degs$symbol %in% up_genes]
down_degs <- degs$symbol[degs$symbol %in% down_genes]

# Add the up and down DEStatus info to starting info
up_v_down_status <- data.frame(up_deg = mix_deg_info$symbol%in% up_degs, down_deg = mix_deg_info$symbol %in% down_degs,
                               log2fc = subj_log2fcs)
all_deg_info <- cbind(mix_deg_info, up_v_down_status)

# Figure out the counts for the respective pathways for up only and down only
  # For bidirectional figure out the counts
  both_status <- all_deg_info$DEstatus
  names(both_status) <- all_deg_info$symbol
  mix_both_counts <- find_contingency_table_counts_given_status(gene_status = both_status, 
                                                                genes_in_go_list = genes_in_go_list)
  
  # For up regulated figure out the counts
  up_status <- all_deg_info$up_deg
  names(up_status) <- all_deg_info$symbol
  mix_up_counts <- find_contingency_table_counts_given_status(gene_status = up_status, genes_in_go_list = genes_in_go_list)
  
  # For down regulated figure out the counts
  down_status <- all_deg_info$down_deg
  names(down_status) <- all_deg_info$symbol
  mix_down_counts <- find_contingency_table_counts_given_status(gene_status = down_status, genes_in_go_list = genes_in_go_list)
  
# Combine all the info and output
subj_direct_info <- list(all_deg_info = all_deg_info, mix_up_counts = mix_up_counts, mix_down_counts = mix_down_counts,
                         mix_both_counts = mix_both_counts)
return(subj_direct_info)
}

find_contingency_table_counts_given_status <- function(gene_status, genes_in_go_list){
#Description: Given a status, it splits the genes into contingency table counts based on the provided status and whether the gene
              #belongs to a given GO term. 
#Input: gene_status (vector of booleans) - provides the positive or negative status for all the genes. Must have entries which
                                          #are named by the genes they represent. 
       #genes_in_go_list (list of vectors of strings) - provides the names of all the genes belonging to each GO term. Each object
                                                        #corresponds to a different GO term. Note that gene names used must match
                                                        #those used for in gene_status names. 
#Output: cont_table_counts (matrix of doubles) - provides a matrix of contingency table counts for each GO term where margins 
                                              #separate genes based on positive status and whether they are present in GO term. 
                                              #Rows = GO terms and 1st column is positive and in GO term, 2nd is negative and in
                                              #GO term, 3rd is positive and outside GO term, and 4th is negative and outside GO 
                                              #term. 
  
# Loop through all the go terms
cont_table_counts <- t(sapply(genes_in_go_list, function(x){
   
  # Find the number of genes in the GO term
  num_genes <- length(gene_status)
  num_genes_in_go <- sum(names(gene_status) %in% x) # Total number of genes which are in the go term
  num_genes_out_go <- num_genes - num_genes_in_go
  num_pos_genes <- sum(gene_status == TRUE)
  num_negative_genes <- num_genes - num_pos_genes
  
  # Find contingency table counts for said GO term
  genes_in_go <- gene_status[names(gene_status) %in% x]
  num_pos_genes_in_go <- sum(genes_in_go == TRUE)                    # Positive in GO term
  num_neg_genes_in_go <- num_genes_in_go - num_pos_genes_in_go    # Negative in GO term
  num_pos_genes_out_go <- num_pos_genes - num_pos_genes_in_go     # Positive Outside GO term
  num_neg_genes_out_go <- num_genes_out_go - num_pos_genes_out_go # Negative Outside GO term
  
  # Combine and name appropriately all the counts and output
  cont_counts <- c(num_pos_genes_in_go, num_neg_genes_in_go, num_pos_genes_out_go, num_neg_genes_out_go)
  names(cont_counts) <- c('altered_inside', 'unaltered_inside', 'altered_outside', 'unaltered_outside')
  return(cont_counts)
}) )
rownames(cont_table_counts) <- names(genes_in_go_list)

return(cont_table_counts)
}

find_genes_in_gos <- function(annofile){
#Description: It provides a list of genes belonging to all the mentioned GO terms. 
#Input: annofile (dataframe) - shows which genes belong to which GO terms. Must have variable called 'path_id' which shows the 
                              #GO terms and the variable 'symbol' which lists the genes belonging to said GO term. 
#Output: genes_in_go_list (list of vectors of strings) - each objects corresponds to a unique GO term and is named for such. Each
                                                        #object contains the genes belonging to said GO term.
  
# If the annotation file contains only 1 unique term
if(length(unique(annofile$path_id)) == 1){
  # Find the genes within the only gene set
  genes_in_one_go <- unique(annofile$symbol)
  
  # PUt the format into a list object containing a vector (so that it works in later pipelines)
  genes_in_go_list <- list(genes_in_one_go = genes_in_one_go)
  names(genes_in_go_list) <- unique(annofile$path_id)
  return(genes_in_go_list)
}
  
# Find the size of each GO term
GO_anno_table <- as.data.table(annofile)
setkey(GO_anno_table, path_id)
genes_in_go_list <- sapply(unique(GO_anno_table$path_id), function(x){
  return(GO_anno_table[.(x), nomatch = 0L]$symbol)
})
names(genes_in_go_list) <- unique(GO_anno_table$path_id)  

# Output the results
return(genes_in_go_list)
}

calc_multi_counts <- function(genes_info, genes_in_gos_list){
#Description: It calculates the contingency table counts for the 2 x 3 contingency table where the margin for the rows are
              #'in_go' and 'out_of_go', and the margins for the columns are 'up DEG', 'down DEG', 'neither up nor down DEG'. 
#Input: genes_info (dataframe) - contains the info for one subject about which genes are differentially expressed and in which
                                #direction. Must have a variable called 'symbol' which shows the gene name, a variable called 
                                #'up_deg' which shows the genes which are up regulated and DEG, and a variable called 'down_deg' 
                                #which shows which genes are down regulated and DEG.  
       #genes_in_gos_list (list of vectors of strings) - lists the genes within each GO term by name. Each object corresponds to
                                                        #a different GO term. 
#Output: multinomial_counts (matrix of integers) - the 2x3 contingency table counts for the subject for each GO term. 
 
# Loop through all the go terms
multinomial_counts <- t(sapply(genes_in_gos_list, function(x){ 
  
  # Find the total number of genes for margins
  num_genes <- nrow(genes_info)
  num_genes_in_go <- length(x)
  num_genes_out_go <- num_genes - num_genes_in_go
  num_pos_genes <- sum(genes_info$up_deg)
  num_neg_genes <- sum(genes_info$down_deg)
  
  # Find genes in GO term
  genes_in_go <- genes_info[genes_info$symbol %in% x, ] # Note it is possible for their to be less genes than in GO term
  
  # Find the individual counts for 1st row of table
  num_pos_in_go <- sum(genes_in_go$up_deg)
  num_neg_in_go <- sum(genes_in_go$down_deg)
  num_neither_in_go <- num_genes_in_go - num_pos_in_go - num_neg_in_go
  
  # Find individual counts for 2nd row of table
  num_pos_out_go <- num_pos_genes - num_pos_in_go
  num_neg_out_go <- num_neg_genes - num_neg_in_go
  num_neither_out_go <- num_genes_out_go - num_pos_out_go - num_neg_out_go 
  
  # Put the counts together and output
  contingency_table_counts <- c(num_pos_in_go, num_neg_in_go, num_neither_in_go, num_pos_out_go, num_neg_out_go, 
                                num_neither_out_go)
  names(contingency_table_counts) <- c('num_pos_in_go', 'num_neg_in_go', 'num_neither_in_go', 'num_pos_out_go', 'num_neg_out_go', 
                                       'num_neither_out_go')
  return(contingency_table_counts)
}) )

# Name and output the counts
return(multinomial_counts)
}