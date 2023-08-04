#' 'compute N-of-1-MixEnrich'
#'
#' @param abs.logFC #abslogFC is calculated as |log2(X1+1)-log2(X2+1)|, where X1 is the expression value of a gene under condition 1, and X2 is the expression values of a gene under condition 2
#' @param case A vector (paired with "baseline") containing gene expression from a sample of interest
#' @param gene.symbol gene symbol is a vector of gene names corresponding to abs.logFC
#' @param GeneSet.annotation GeneSet.list is a list of genesets. Each element of this geneset contains the gene names of the genes belonging to the gene set.
#' @param FDR.method method for multiplicity correction
#' @param alternative indicates the alternative hypothesis in Fisher's exact test.
#' 
#' @return A list contains two elements. The first element is the pathways and their corresponding p.values. The second element is all the genes and their corresponding DEG status
#' If either input vectors have no variance, the MD distance is chosen to the simple difference.

mixEnrich = function (abs.logFC, gene.symbol,GeneSet.annotation, FDR.method = 'BY', FDR=T,alternative = c('two.sided','greater','less'),
                      abs.log2FC_Cutoff =0,posterProb){ 

set.seed(1234)
alternative = match.arg(alternative)
abs.logFC.tmp = abs.logFC[abs.logFC != 0 ]
mix1 = try(normalmixEM(abs.logFC.tmp,lambda=.5,mu=c(.5,3),sigma=c(1,1)) )
p.dereg = mix1$posterior[,which.max(mix1$mu)]
  ##find the genes that we  believe they are deregulated
ind.DE = logical(length= length(abs.logFC)) 
  #Above line: Creating a logical vector of same length as abs.logFC

ind.DE[abs.logFC != 0 ] = (p.dereg>posterProb)&(abs.logFC.tmp > abs.log2FC_Cutoff)
  #Above line: Set all indices (where abs.logFC !=0) as 1 if posterior > .5 for being dysregulated
  #So, above line is list of 1 or 0 if gene is dysregulated and it is list for all genes
  
posterior_probs <- rep(0, length(abs.logFC))
posterior_probs[abs.logFC != 0 ] = p.dereg
DE.df = data.frame(symbol = gene.symbol,abslogFC = abs.logFC, posterior.probability = posterior_probs,
                   DEstatus = as.integer(ind.DE))

GeneSet.annotation = GeneSet.annotation[GeneSet.annotation$symbol %in% gene.symbol,]
GO.sym.DE = merge(GeneSet.annotation,DE.df,by = 'symbol')
## The total number of genes in each GO category (is t1)
t1 = table(GO.sym.DE$path_id) 
## the number of DE genes in each GO (is t2)
list2 <- lapply(unique(GO.sym.DE$path_id), function(x) { 
  gene_set <- GO.sym.DE[GO.sym.DE$path_id == x, ]
  return(data.frame(GO_ID = gene_set$path_id[1], DE = sum(gene_set$DEstatus)) )
})
df.GO2 <- do.call( 'rbind', list2)
df.GO1 = as.data.frame(t1)
names(df.GO1) = c('GO_ID', 'total')
## GO.count is a data frame which contains GO ids and the total number of genes belongs to this GO term and the DEGs in this GO term
  GO.count = merge(df.GO1, df.GO2, by = 'GO_ID')
    
  ## background gene DE status
  bg.DE  = c(sum(ind.DE), length(ind.DE) - sum(ind.DE))
    #Above line: bg.DE = (Num DEG over all genes, Num not DEG over all genes)

    #Create a matrix to contain the Odds ratios
     oddsRatio<-matrix(0,nrow = dim(GO.count)[1],ncol = 1)
     
#Create an array to hold all the contingency tables
 contTables<-array(0,c(2,2,dim = dim(GO.count)[1])) 
 
    p.val.fisher = numeric(dim(GO.count)[1])      # a vecter to store the p.values of FET for each GO
    for (i in 1:dim(GO.count)[1]){ # calculate the p.value of fisher exact tests for every gene set  
       single.GO = rbind(GO.count[i,]$DE,GO.count[i,]$total-GO.count[i,]$DE)
         #Above line: single.GO = (Num DEG in pathway , Num not DEG in pathway)

        m = cbind(single.GO, bg.DE - single.GO)
         #So m is 2 by 2 matrix: (Num DEG in pathway, Num DEG not in pathway) over (Num not DEG in pathway, Num not DEG not in pathway)

        m = t(m)
         #Above line: Transposes m, so now is (Num DEG in pathway, Num not DEG in pathway) 
          #over (Num DEG not in pathway, Num not DEG not in pathway)
 
        contTables[,,i]<-m #Load the contingency table into the array
        
        p.val.fisher[i]=fisher.test(m, alternative = alternative)$p.value
         #Above line: Obtains the p-value for the fisher test

        #Obtain the estimate of the odds ratio 
        oddsRatio[i,1]<-fisher.test(m,alternative = alternative)$estimate
} 
if(FDR){
 #Find the lnOR and variance for the lnOR
  ai = contTables[1,1,]
  bi = contTables[1,2,]
  ci = contTables[2,1,]
  di = contTables[2,2,]
  lnOR<-log( ((ai+.5 )*(di+.5) )/((bi+.5)*(ci+.5)) )
  varlnOR<- (1/(ai+.5)) + (1/(bi+.5))+(1/(ci+.5))+(1/(di+.5))
  
 #Adjust the fisher exact test p-value 
    p.val.fisher.adj = p.adjust(p.val.fisher,method=FDR.method)
    fisher.res = data.frame(GO_ID = GO.count$GO_ID, p.value.adjusted= p.val.fisher.adj, p.value = p.val.fisher, Odds_Ratio = oddsRatio,
                            ai = contTables[1,1,],bi = contTables[1,2,],ci = contTables[2,1,], di = contTables[2,2,],
                            lnOR = lnOR, varlnOR = varlnOR)
    return (c(fisher.res,DE.df))
}else{
    fisher.res = data.frame(GO_ID = GO.count$GO_ID, p.value = p.val.fisher)
return (fisher.res)
}
}

library(mixtools)
## define a function to calculate absolute log Fold change
absLogFC = function(x1,x2){
    abs(log2(x1+1) - log2(x2+1))
}
