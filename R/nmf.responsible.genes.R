#' A Function to find responsible genes for variation within a variable
#'
#' This function allows you to find the likely genes responsible for affecting a variable, typically the presence or absence of a disease
#' @param H.matrix The H.matrix of an NMF factorisation of a gene expression dataset
#' @param W.matrix The W.matrix of an NMF factorisation of a gene expression dataset
#' @param variable The discrete variable or disease to be investigated
#' @param number.genes The number of likely genes to find
#' @param gene.names A character/string vector containing the list of gene names
#' @param gene.expressions The initial gene expression matrix, with genes labelled by row names
#' @keywords cats
#' @export
#' @examples
#' nmf.responsible.genes()









nmf.responsible.genes <- function(H.matrix,W.matrix,variable,number.genes,gene.names,gene.expressions,control.age = T,age = NA){
  if (control.age == T){
  responsible.meta.gene = which.max(abs(cor(as.numeric(as.factor(variable)),t(H.matrix))) - 
                                          abs(cor(as.numeric(as.factor(age)),t(H.matrix))))
  } else {
  responsible.meta.gene = which.max(abs(cor(as.numeric(as.factor(variable)),t(H.matrix))))
}
  top.gene_IDs <- (names(sort(W.matrix[,responsible.meta.gene],decreasing=TRUE)[1:number.genes]))
  top.gene.names <- gene.names[pmatch(top.gene_IDs,row.names(gene.expressions))]
  return(top.gene.names)
}