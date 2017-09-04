gene.expression.heatmap <- function(gene.expression.matrix, number.genes, log.scaled,row.scaled){
  
  
  if (log.scaled == TRUE){
    gene.expression.matrix <- log1p(gene.expression.matrix)
  }
  if (row.scaled == TRUE){
    gene.expression.matrix <- gene.expression.matrix/(rowMax(gene.expression.matrix)+1e-10)
  }
  varying.row.index <- RowCV(gene.expression.matrix)
  gene.expression.matrix.sample <- 
    gene.expression.matrix[order(varying.row.index,decreasing=T)[seq(1,nrow(gene.expression.matrix),nrow(gene.expression.matrix)/number.genes)],]
  heatmap.2(gene.expression.matrix.sample,
            col=viridis::viridis(200),
            dendrogram='none',
            Rowv=FALSE,
            Colv=FALSE,
            trace='none',
            density.info="none",
            labCol=FALSE,
            labRow=FALSE)
  
}