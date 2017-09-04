#' A Function to find responsible genes for variation within a variable
#'
#' This function allows you to find the likely genes responsible for affecting a variable, typically the presence or absence of a disease
#' @param k.range Suggest a vector of k values to explore, with a default
#' @param zero.threshold what proportion of zeros should a gene have before being left out of the analysis
#' @param variable The discrete variable or disease to be investigated
#' @param number.genes The number of likely genes to find
#' @param gene.names A character/string vector containing the list of gene names
#' @param gene.expressions The initial gene expression matrix, with genes labelled by row names
#' @keywords cats
#' @export
#' @examples
#' nmf.responsible.genes()

nmf.gene.choose.k <- function(k.range = seq(2,ncol(gene_expressions)/10,2), 
                              gene_expressions, 
                              zero.threshold = 0.1, 
                              log.scaled = T,
                              loss.type = 'mse',
                              b=c(0,0,0),
                                na = FALSE,
                                sample = FALSE,
                                number.genes,a=c(0,0,0)) {
  require(NNLM)
  zeros_geneRpkm <- gene_expressions==0
  
  # how many zero entries per gene; look to reduce features
  zeros_geneRpkm_per_gene <- rowSums(zeros_geneRpkm)
  
  # how many genes have at least one zero entry
  sum(zeros_geneRpkm_per_gene>0)
  
  # about 45,000 have at least one zero entry, so can remove 45,000 entries, leaving approx 18,000
  # now create IDs of genes to be eliminated
  
  zero_gene_index <- which(zeros_geneRpkm_per_gene>zero.threshold*dim(gene_expressions)[2])
  
  #prune original gene data
  print(paste(length(zero_gene_index), "dimensions removed"))
  if (length(zero_gene_index) > 0){
    pruned_expression_data <- gene_expressions[-zero_gene_index,]
  } else {
    print("activated")
    pruned_expression_data <- gene_expressions
  }
  if (na == TRUE){
    pruned_expression_data[pruned_expression_data ==0] <- NA
  }
  if (log.scaled == TRUE){
    print("yo")
    pruned_expression_data <- log1p(pruned_expression_data)
  }
  if (sample == TRUE){
    varying.row.index <- RowCV(pruned_expression_data)
    pruned_expression_data <- pruned_expression_data[order(varying.row.index,decreasing=T)[1:number.genes],]
  }
  error.df <- data.frame(matrix(nrow = length(k.range), ncol = 4))
  colnames(error.df) <- c('k','mse loss','mkl loss','target loss')
  j = 1
  for (i in k.range){
    nmf_data <- nnmf(pruned_expression_data,n.threads=0,k=i,loss = loss.type, beta = b, alpha = a)
    error.df[j,] <- c(i,nmf_data$mse,nmf_data$mkl,nmf_data$target.loss)
    j = j+1
  }
  return(error.df)

}