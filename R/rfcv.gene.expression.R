rfcv.gene.expression <- function(folds, step,
                                 sample = FALSE,
                                 number.genes,
                                 gene_expressions,
                                 log.scaled,
                                 zero.threshold,
                                 variable) {
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
  if (log.scaled == TRUE){
    print("yo")
    pruned_expression_data <- log1p(pruned_expression_data)
  }
  if (sample == TRUE){
    varying.row.index <- RowCV(pruned_expression_data)
    pruned_expression_data <- pruned_expression_data[order(varying.row.index,decreasing=T)[1:number.genes],]
  }
  rfcv_data <- rfcv(train=t(pruned_expression_data),trainy=as.factor(variable),step = step, cv.fold = folds)
}