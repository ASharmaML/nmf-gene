nmf.boxplot <- function(H.matrix,variable,pd,plot.layout){
  class.data <- list()
  kcount = 1
  par(mfrow=plot.layout)
  for (i in 1:dim(H.matrix)[1]){
    
    scaled_H <- H.matrix[kcount,]/max(H.matrix[kcount,]+0.00000000001)
    class.data[[1]] <- scaled_H[variable==variable[1]]
    class.data[[2]] <- scaled_H[variable!=variable[1]]
    kcount <- kcount + 1
    
    boxplot(class.data, ylab = "Values",xlab = paste("Dimension",i), notch = FALSE, col = c("red","dark green"), xaxt='n')
    legend(x=1,y=1, ncol=1,  cex=1, c(variable[1],variable[variable!=variable[1]][1]),text.width=5,lty=c(1,1),lwd=c(6,6),col=c("red","dark green"))
  }
  
  #axis(side=1, at=seq(1,96,2),labels=c(1:48))
  #legend(x=1,y=600, ncol=1,  cex=1, c("Class 1","Class 2"),text.width=5,lty=c(1,1),lwd=c(6,6),col=c("red","dark green"))
  #boxplot(data.df2[,2:29], ylab = "Values", notch=TRUE)
  title(main="All Box Plots", line = 2,xlab="Dimension Number",  ylab="Values", outer=TRUE)
}