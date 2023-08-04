runMixEnrich <- function( theAbsLogFCs, subNames, myPosterProb, myAlternative, GO_anno, absLogFC_Cutoff = 0, num_core = 4 ){
#Description: This function runs an inputted version of mixEnrich on all the pairs of data and outputs the results. The results
              #show both gene set level results and transcript level results. 
#Inputs: theAbsLogFCs (matrix of doubles) - the absolute log2FCs of the tumor over normal paired expression values
        #subNames (character vector) - the names you want to assign the list objects in the output. Ideally these would match
                                      #up or be equal to the column names of 'theAbsLogFCs' input. These were meant to be 
                                      #subject names. 
        #myPosterProb (double) - the posterior probability cutoff for determining altered transcript status for MixEnrich
        #alternative (string) - either 'two.sided','greater',or 'less' it is the alternative hypothesis you wish the fisher exact test in 
                               #MixEnrich to test
        #GO_anno (dataframe) - shows which genes/transcripts belong to which gene sets. Needs to have two variables. The 1st 
                              #variable must be called 'path_id' and must show the gene set that a gene/transcript belongs to. 
                              #The 2nd variable must be called 'symbol' and shows what the name of the gene/transcript is. 
        #absLogFC_Cutoff (double) - the cutoff for the abs(log2FC) which is used for determining altered transcript status for
                                    #MixEnrich. 
        #num_core (integer) - the number of cores you want to call for parallel running of MixEnrich. Needs to be a positive
                              #integer which is less than or equal to the number of cores available on your computer. 
#Outputs: myOutput (list) - first object is a list of containing all the MixEnrich results for each individual and 2nd object is 
                           #matrix of altered transcript statuses where each column corresponds to results for each individual
  
#Create the desired cluster
cl <- makeCluster(num_core)
clusterExport(cl, c('theAbsLogFCs', 'myPosterProb', 'myAlternative', 'GO_anno', 'mixEnrich', 'absLogFC_Cutoff'),
              envir = environment())
clusterEvalQ(cl, library(mixtools))
  
#Loop through all the patients 
mixRes <- parApply(cl, theAbsLogFCs, 2, myFun <- function(i)
{
  #Run the mixEnrich on the pair (Note that the alt is two-sided)
  tempRes <- mixEnrich(abs.logFC = i, gene.symbol = names(i), 
                     GeneSet.annotation = GO_anno, FDR.method = 'BH',
                     alternative = myAlternative, posterProb = myPosterProb, abs.log2FC_Cutoff = absLogFC_Cutoff)
  
  #Save the data in dataframe format
  pathwayInfo <- data.frame(GO_ID = tempRes$GO_ID, p.value.adjusted = tempRes$p.value.adjusted, p.value = tempRes$p.value,
                          Odds_Ratio = tempRes$Odds_Ratio, DEG_in_Pathway = tempRes$ai,
                          not_DEG_in_Pathway= tempRes$bi, 
                          DEG_not_in_Pathway = tempRes$ci, 
                          not_DEG_not_in_Pathway = tempRes$di, lnOR = tempRes$lnOR, varlnOR = tempRes$varlnOR)
    
  #Save the gene data into the DEGs matrix 
  DEGInfo <- data.frame(symbol = tempRes$symbol, absLogFC = tempRes$abslogFC,
                      posterior.probability = tempRes$posterior.probability,
                      DEstatus = tempRes$DEstatus)
  rownames(DEGInfo) <- names(theAbsLogFCs)
  
  #Return the combined pathway and DEG results
  myOut <- vector('list', 2)
  myOut[[1]] <- pathwayInfo
  myOut[[2]] <- DEGInfo
  return(myOut)
})  
#End of loop through all the patients
  
#Close the cluster you created to return the memory
stopCluster(cl)
  
#Combine all of the pathway info across patients into a single list and all the DEG info into a single list dim(Data)[2]/2
nOf1Results <- lapply(1:(ncol(theAbsLogFCs)), function(i){ mixRes[[i]][[1]] }) # Combine pathway info
names(nOf1Results) <- subNames  # Label the results by the individual they belong to

DEGs <- lapply(1:(ncol(theAbsLogFCs)), function(i){ mixRes[[i]][[2]]}) #Combine transcript info into a matrix
names(DEGs) <- subNames #Name the columns by their corresponding subject ids
  
#Combine the results and output them
myOutput <- list(nOf1Results = nOf1Results, 
                 DEGs = DEGs)
return(myOutput)
}
#End of function for running MixEnrich

obtain_mixenrich_counts <- function(log2fcs, subj_names, poster_prob_cutoff, the_alternative, anno_file, 
                                    abs_log2fc_cutoff = 0, num_core = 4){
#Description: It runs the MixEnrich pipeline and then finds the mixenrich contingency table counts for each of the directions:
              #ie both, up, down. 
#Input: log2fcs (matrix of doubles) - the log2FCs corresponding to the within-subject paired samples. Needs to have columns
                                      #named after the subjects and said names should match subj_names input. Rows should be 
                                      #named by GO term. (or transcription factor, whichever term is under path_id in anno_file). 
       #subj_names (vector of strings) - the names for the subjects. Must match the column names of log2FCs. 
       #poster_prob_cutoff (double) - the posterior probability cutoff for determining DEGs. Must be between 0 and 1. 
       #the_alternative (string) - the type of alternative hypothesis you want for MixEnrich. Possible options are 'two.sided', 
                                  #'greater', and 'less'. 
       #anno_file (dataframe) - contains the GO terms (or transcription factors) and what genes belong to them. Must have a 
                                #variable called 'path_id' which lists the GO terms (or TFs) and must have a variable called 
                                #'symbol' which provides the genes which belong to the GO term (or TF). 
       #abs_log2fc_cutoff (double) - the |log2FC| cutoff for determining DEG status in MixEnrich. Should be positive. 
       #num_core (integer) - the number of cores you want to call for parallel running of MixEnrich. Needs to be a positive
                            #integer which is less than or equal to the number of cores available on your computer. 
#Output: count_output (list of matrices and dataframe) - Contains the mixenrich counts and the MixEnrich output.
  
# Run MixEnrich to obtain gene set level and gene level results
tempResults <- runMixEnrich(theAbsLogFCs = abs(log2fcs), subNames = subj_names, myPosterProb = poster_prob_cutoff,
                            myAlternative = the_alternative, GO_anno = anno_file, absLogFC_Cutoff = abs_log2fc_cutoff, 
                            num_core = num_core)
nOf1Results <- tempResults[[1]] #Enter the Pathway level info from the results
DEGs <- tempResults[[2]] #Enter the DEG info from the results 
  
# Add upregulation vs downregulation for the MixEnrich info 
  # Find the list of genes within each GO term
  genes_in_gos_list <- find_genes_in_gos(annofile = anno_file)

  # For each subject figure out the contingency table counts
  subj_de_dir_info <- lapply(1:length(DEGs), function(i){
  return(add_side_regulation_info(mix_deg_info = DEGs[[i]], 
                                  subj_log2fcs = log2fcs[, colnames(log2fcs) == names(DEGs)[i]],
                                  genes_in_go_list = genes_in_gos_list))
  })

# Combine the info types into one
complete_mix_deg_info <- lapply(subj_de_dir_info, function(x){ return(x$all_deg_info) })
mixenrich_counts_up <- lapply(subj_de_dir_info, function(x){ return(x$mix_up_counts) })
mixenrich_counts_down <- lapply(subj_de_dir_info, function(x){ return(x$mix_down_counts) })
mixenrich_counts_both <- lapply(subj_de_dir_info, function(x){ return(x$mix_both_counts) })
names(complete_mix_deg_info) <- names(mixenrich_counts_up) <- names(mixenrich_counts_down) <- names(DEGs)
names(mixenrich_counts_both) <- names(DEGs)

# Convert the MixEnrich nof1 Results into matrices of counts
mixenrich_counts <- lapply(nOf1Results, function(x){
  the_counts <- x[, c(5:8)] 
  colnames(the_counts) <- c('altered_inside', 'unaltered_inside', 'altered_outside', 'unaltered_outside')
  rownames(the_counts) <- x$GO_ID
  return(the_counts)
})

#Combine all of the results and output
counts_output <- list(mixenrich_counts = mixenrich_counts, mixenrich_counts_up = mixenrich_counts_up,
                      mixenrich_counts_down = mixenrich_counts_down, nOf1Results = nOf1Results, 
                      complete_mix_deg_info = complete_mix_deg_info)
return(counts_output)
}

runMixEnrich_slow <- function( theAbsLogFCs, subNames, myPosterProb, myAlternative, GO_anno, absLogFC_Cutoff = 0){
#Description: This function runs an inputted version of mixEnrich on all the pairs of data and outputs the results without using
              #parallel code. 
#Inputs: theAbsLogFCs (matrix of doubles) - the absolute log2FCs of the tumor over normal paired expression values
        #subNames (character vector) - the character list of the comparison names corresponding to the individuals
        #myPosterProb (double) - the posterior probability cutoff for determining DEG status for MixEnrich
        #alternative (string) - either 'two.sided','greater',or 'less' it is the alternative hypothesis you wish the fisher exact test in 
                               #MixEnrich to test
        #GO_anno (dataframe) - the file which shows which genes/transcripts belong to which gene sets. Needs to have two 
                              #variables. The 1st variable must be called 'path_id' and must show the gene set that a 
                              #gene/transcript belongs to. The 2nd variable must be called 'symbol' and shows what the name of
                              #the gene/transcript is. 
        #absLogFC_Cutoff (double) - the cutoff for the abs(log2FC) which is used for determining DEG status for MixEnrich. 
#Outputs: myOutput (list) - first object is a list of containing all the MixEnrich results for each individual and 2nd object is 
                            #matrix of DEG statuses where each column corresponds to results for each individual
  

#Loop through all the patients 
mixRes <- apply(theAbsLogFCs, 2, function(i){
  #Run the mixEnrich on the pair (Note that the alt is two-sided)
  tempRes<-mixEnrich(abs.logFC = i, gene.symbol = names(i), 
                     GeneSet.annotation = GO_anno, FDR.method = 'BH',
                     alternative=myAlternative,posterProb = myPosterProb,abs.log2FC_Cutoff =absLogFC_Cutoff)
    
  #Save the data in dataframe format
  pathwayInfo<-data.frame(GO_ID=tempRes$GO_ID,p.value.adjusted=tempRes$p.value.adjusted,p.value=tempRes$p.value,
                          Odds_Ratio = tempRes$Odds_Ratio,DEG_in_Pathway = tempRes$ai,
                          not_DEG_in_Pathway= tempRes$bi, 
                          DEG_not_in_Pathway = tempRes$ci, 
                          not_DEG_not_in_Pathway= tempRes$di,lnOR = tempRes$lnOR, varlnOR = tempRes$varlnOR)
    
  #Save the gene data into the DEGs matrix 
  DEGInfo<-data.frame(symbol = tempRes$symbol, absLogFC = tempRes$abslogFC,
                      posterior.probability = tempRes$posterior.probability,
                      DEstatus = tempRes$DEstatus)
  rownames(DEGInfo)<-names(theAbsLogFCs)
  
  #Return the combined pathway and DEG results
  myOut<-vector('list',2)
  myOut[[1]]<-pathwayInfo
  myOut[[2]]<-DEGInfo
  return(myOut)
})  
#End of loop through all the patients
  
#Combine all of the pathway info across patients into a single list and all the DEG info into a single list dim(Data)[2]/2
nOf1Results <- lapply(1:(dim(theAbsLogFCs)[2]), function(i){ mixRes[[i]][[1]] }) #Combine pathway info
names(nOf1Results) <- subNames #Label the results by the individual they belong to
DEGs <- lapply(1:(dim(theAbsLogFCs)[2]), function(i){ mixRes[[i]][[2]]}) #Combine DEG info into a matrix
names(DEGs) <- subNames #Name the columns by their corresponding subject ids
  
#Combine the results and output them
myOutput<-list(nOf1Results,DEGs)
return(myOutput)
}
#End of function for running MixEnrich

