#' A Function to find responsible genes for variation within a variable
#'
#' This function allows you to find the likely genes responsible for affecting a variable, typically the presence or absence of a disease
#' @param H.matrix The H.matrix of an NMF factorisation of a gene expression dataset
#' @param variable A vector containing values for the discrete variable or disease to be investigated
#' @param step The skip value to determine how many patients/subjects to display
#' @param variable.name The name of the variable being investigated
#' @keywords cats
#' @export
#' @examples
#' nmf.responsible.genes()

nmf.heatmap <- function(H.matrix,variable,step=1,variable.name = "unspecified variable"){
  if (is.numeric(variable)){
    colsep.vec <- F
    labvec.values <- 
      variable[order(as.vector(variable))[seq(1,length(variable),(10*step))]]
    
    labvec <- rep(NA,length(variable)/step)
    
    labvec[seq(1,length(variable)/step,10)] <- round(labvec.values,digits=1)[1:length(seq(1,length(variable)/step,10))]
  } else {
    colsep.vec <- c(table(variable)[[1]]/step)
    labvec <- rep(NA,length(variable)/step)
    sum.factors <- 0
    for (i in 1:length(levels(as.factor(variable)))){
      labvec[(table(variable)[i]/2+sum.factors)/step] <- levels(as.factor(variable))[i]
      sum.factors = sum.factors + table(variable)[i]
    }
   
  }

  

  par(mar=c(1,1,1,1))
  H <- (H.matrix - rowMin(H.matrix))/(rowMax(H.matrix-rowMin(H.matrix))+1e-10)
  heatmap.2(H[,order(as.vector(variable),
                     max.col(t(H)))[seq(1,dim(H)[2],step,)]],
            col=viridis::viridis(200),dendrogram='none',
            Rowv=FALSE,Colv=FALSE,trace='none',
            colsep=colsep.vec,
            sepwidth=c(1,0.05),
            labCol = labvec,
            xlab = paste("Subjects by",variable.name),
            ylab = "Meta-Gene Expressions",
            main = "Meta-Gene Expression of Patients",
            density.info = "none",
            keysize = 1)
  par(mar=c(2,2,2,2))

}