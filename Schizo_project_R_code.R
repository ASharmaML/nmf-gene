# Project code for schrizophenia dataset 
source("https://bioconductor.org/biocLite.R")
biocLite()
library("NNLM")
biocLite("affy")
library("affy")
library("NMF")
library("corrplot")
library("plyr")
library("tsne")
library("gplots")
library("randomForest")
biocLite("limma")
######### load data in
library("viridis")
library("raster")
library("nmfgene")
library("ggplot2")
library("reshape2")
library("rbenchmark")
biocLite("multtest")
library("multtest")
library("devtools")
library("autoencoder")
setwd("C:/Users/aseem/OneDrive/dropbox/UCL_ML/Term_3/Project/Data/Schizophrenia")


### sample leukemia data
leukemia_big <- read.csv("http://web.stanford.edu/~hastie/CASI_files/DATA/leukemia_big.csv")



load("Schizo_project_data.RData")

load("geneLevel_LIBD_qSVsAndMod.rda") ## first dataset
geneRpkm.1 <- geneRpkm
geneMap.1 <- geneMap
mod.1 <- mod
qSVs.1 <- qSVs
pd.1 <- pd





load("geneLevel_CMC_qSVsAndMod.rda") ## second data set
geneRpkm.2 <- geneRpkm
geneMap.2 <- geneMap
mod.2 <- mod
qSVs.2 <- qSVs
pd.2 <- pd

#### remove rows with all 0s



######################

######## Investigate data/Explore ##############

# variable explanations
### geneMap contains information about each gene
# geneRpkm contains all the schizo patient data
# pd contains additional information about the patients
# qSVs contains the 12 principal components



######## investigate zero entries

# create id matrix of zeros

zeros_geneRpkm <- geneRpkm.1==0

# how many zero entries per gene; look to reduce features
zeros_geneRpkm_per_gene <- rowSums(zeros_geneRpkm)

# how many genes have at least one zero entry
sum(zeros_geneRpkm_per_gene>0)

# about 45,000 have at least one zero entry, so can remove 45,000 entries, leaving approx 18,000
# now create IDs of genes to be eliminated

zero_gene_index <- which(zeros_geneRpkm_per_gene>0)

#prune original gene data

pruned_expression_data.1 <- geneRpkm.1[-zero_gene_index,]
pruned_expression_data.1 <- log1p(pruned_expression_data)


# now for the fun nmf stuff
# first log scale the data ?? Already log scaled??

# classify according to gender etc
rfcv_nnmf_schizo <- rfcv( trainx=t(nmf_20_df[1:20,]),trainy=as.factor(pd$Dx),ntree=500,cv.fold=5,step=0.95)

rfcv_nnmf_schizo_12 <- rfcv( trainx=t(nmf_data_10$H[1:10,]),trainy=as.factor(pd$Dx),ntree=500,cv.fold=5,step=0.95)


rfcv_pca_schizo <- rfcv( trainx=qSVs[,1:12],trainy=as.factor(as.logical(nmf_20_df[21,])),ntree=500,cv.fold=5,step=0.95)

rfcv_nnmf_gender <- rfcv( trainx=t(nmf_20_df[1:20,]),trainy=as.factor(pd$Sex),ntree=500,cv.fold=5,step=0.95)
rfcv_pca_gender <- rfcv( trainx=qSVs[,1:12],trainy=as.factor(pd$Sex),ntree=500,cv.fold=5,step=0.95)

rfcv_nnmf_race <- rfcv( trainx=t(nmf_20_df[1:20,]),trainy=as.factor(pd$Race),ntree=500,cv.fold=5,step=0.95)
rfcv_pca_race <- rfcv( trainx=qSVs[,1:12],trainy=as.factor(pd$Race),ntree=500,cv.fold=5,step=0.95)

rfcv_nnmf_smoke <- rfcv( trainx=t(nmf_20_df[1:20,]),trainy=as.factor(pd$SmokingEither),ntree=500,cv.fold=5,step=0.95)
rfcv_pca_smoke <- rfcv( trainx=qSVs[,1:12],trainy=as.factor(pd$SmokingEither),ntree=500,cv.fold=5,step=0.95)
## Independent Component Analysis

# plots
plot(rfcv_nnmf_schizo$error.cv,rfcv_pca_schizo$error.cv)

# ENSG0000210082


# order original data by race
nmf_data_2_race <- nnmf(pruned_expression_data,n.threads=0,k=20)

heatmap.2(nmf_data_20$W[sample(45101,20),],col=redgreen(75),dendrogram='none',Rowv=FALSE, Colv=FALSE,trace='none')
heatmap.2(nmf_data_20$H[,order(as.vector(pd$Dx))][1:20,seq(1,351,5)],col=redgreen(75),dendrogram='none',Rowv=FALSE, 
          Colv=FALSE,trace='none')
heatmap.2(nmf_data_20$H[,order(as.vector(pd$Race))][1:20,seq(1,351,5)],col=redgreen(75),dendrogram='none',Rowv=FALSE, 
          Colv=FALSE,trace='none')
heatmap.2(nmf_data_20$H[,order(as.vector(pd$Sex))][1:20,seq(1,351,5)],col=redgreen(75),dendrogram='none',Rowv=FALSE, 
          Colv=FALSE,trace='none')
heatmap.2(nmf_data_2_race$H[1:2,seq(1,351,2)],col=redgreen(75),dendrogram='none',Rowv=FALSE, Colv=FALSE,trace='none')


par(mfrow=c(4,5))
for (i in 1:20){
  plot(seq(1,351,1),nmf_data_20$H[i,order(as.vector(pd.1$Sex))])
}

par(mfrow=c(2,3))
for (i in 1:5){
  plot(seq(1,351,1),coef(basic_nmf)[i,order(as.vector(pd.1$Race))])
}


#######################

## Using the two datasets together

# Compare them

summary(rowMax(geneRpkm.2))
summary(rowMax(geneRpkm.1))

geneRpkm.joint <- cbind(geneRpkm.1,geneRpkm.2)
zeros_geneRpkm.joint <- geneRpkm.joint==0

# how many zero entries per gene; look to reduce features
zeros_geneRpkm.joint_per_gene.1 <- rowSums(zeros_geneRpkm.joint[,1:351])
zeros_geneRpkm.joint_per_gene.2 <- rowSums(zeros_geneRpkm.joint[,352:682])

# how many genes have at least one zero entry
sum(zeros_geneRpkm.joint_per_gene>0)

# about 45,000 have at least one zero entry, so can remove 45,000 entries, leaving approx 18,000
# now create IDs of genes to be eliminated

zero_gene_index <- which(zeros_geneRpkm.joint_per_gene.1>0 | zeros_geneRpkm.joint_per_gene.2>0)

#prune original gene data

pruned_joint_data <- geneRpkm.joint[-zero_gene_index,]
scaled_joint_data <- log1p(pruned_joint_data)


nmf_joint_20 <- nnmf(scaled_joint_data,n.threads=0,k=20)


# heatmap

set.seed(10)

heatmap.2(nmf_joint_20$W[sample(dim(nmf_joint_20$W)[1],50),],col=redgreen(75),dendrogram='none',
          Rowv=FALSE, Colv=FALSE,trace='none')
heatmap.2(nmf_joint_20$H[,1:20],col=redgreen(75),dendrogram='none',Rowv=FALSE, Colv=FALSE,trace='none')

pd.joint.Dx <- c(pd.1$Dx,pd.2$Dx)
pd.joint.Dx[pd.joint.Dx=="SCZ"] <- "Schizo"

pd.joint.Gender <- c(pd.1$Sex,pd.2$Gender)
pd.joint.Gender[pd.joint.Gender=="Male"] <- "M"
pd.joint.Gender[pd.joint.Gender=="Female"] <- "F"

pd.joint.Race <- c(pd.1$Race,pd.2$Ethnicity)
pd.joint.Race[pd.joint.Race=="African-American"] <- "AA"
pd.joint.Race[pd.joint.Race=="Caucasian"] <- "CAUC"

par(mfrow=c(4,5))
for (i in 1:20){
  plot(seq(1,dim(nmf_joint_20$H)[2],1),nmf_joint_20$H[i,order(as.vector(pd.joint.Race))])
}


# box plots



nmf_joint_df <- as.data.frame(nmf_joint_20$H)

rfcv_nnmf_gender_joint <- rfcv( trainx=t(nmf_joint_df[1:20,]),trainy=as.factor(pd.joint.Gender),ntree=500,cv.fold=5,step=0.95)

rfcv_nnmf_schizo_joint <- rfcv( trainx=t(nmf_joint_df[1:20,]),trainy=as.factor(pd.joint.Dx),ntree=500,cv.fold=5,step=0.95)

rfcv_nnmf_race_joint <- rfcv( trainx=t(nmf_joint_df[1:20,]),trainy=as.factor(pd.joint.Race),ntree=500,cv.fold=5,step=0.95)


rfcv_pca_gender_joint <- rfcv( trainx=qSVs[,1:12],trainy=as.factor(pd$Sex),ntree=500,cv.fold=5,step=0.95)


# boxplots

class.data.list <- list()
kcount <- 1
for (i in seq(1,39,2))
{
  scaled_H <- nmf_joint_20$H[kcount,]/max(nmf_joint_20$H[kcount,])
  class.data.list[[i]] <- scaled_H[pd.joint.Race=="AA"]
  class.data.list[[i+1]] <- scaled_H[pd.joint.Race!="AA"]
  kcount <- kcount + 1
}
par(mfrow=c(1,1))
boxplot(class.data.list, ylab = "Values",xlab = "Dimension Number", notch = TRUE, col = c("red","dark green"), xaxt='n')
axis(side=1, at=seq(1,96,2),labels=c(1:48))
legend(x=1,y=600, ncol=1,  cex=1, c("Class 1","Class 2"),text.width=5,lty=c(1,1),lwd=c(6,6),col=c("red","dark green"))
#boxplot(data.df2[,2:29], ylab = "Values", notch=TRUE)
title(main="All Box Plots", line = 2,xlab="Dimension Number",  ylab="Values", outer=TRUE)


### compare the 2 different datasets
set.seed(10)
sample.index <- sample(dim(geneRpkm)[1],6)
summary(t(geneRpkm.1[sample.index,]))
summary(t(geneRpkm.2[sample.index,]))



### generalised NMF function for gene expression datasets

RowCV <- function(x) {
  sqrt((rowSums(((x - rowMeans(x))^2)/(dim(x)[2] - 1)))/(rowMeans(x)+1e-10) * 100)
}

nmf.gene.expression <- function(n.clusters = 20, gene_expressions, zero.threshold = 0.1, log.scaled,loss.type = 'mse',b,
                                na = FALSE,
                                sample = FALSE,
                                number.genes,a) {
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
  print(min(pruned_expression_data))
  nmf_data <- nnmf(pruned_expression_data,n.threads=0,k=n.clusters,loss = loss.type, beta = b, alpha = a)
  
}

nmf_15.1 <- nmf.gene.expression(15,geneRpkm.1,zero.threshold = 0.2,log.scaled = TRUE)
par(mfrow=c(3,5))
for (i in 1:15){
  plot(seq(1,351,1),nmf_15.1$H[i,order(as.vector(pd.1$Sex))])
}

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


nmf.boxplot(nmf_15.1$H,pd=pd.1,variable=pd.1$Dx)

nmf.15.1.rfcv.Dx <- rfcv( trainx=t(nmf_15.1$H[1:15,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=5,step=0.95)

nmf_15.2 <- nmf.gene.expression(15,geneRpkm.2,zero.threshold = 0.2, log.scaled =TRUE)
nmf.15.2.rfcv.Dx <- rfcv(trainx=t(nmf_15.2$H[1:15,]),trainy=as.factor(pd.2$Dx),ntree=500,cv.fold=5,step=0.95)
nmf.boxplot(nmf_15.2$H,pd=pd.2,variable=pd.2$Gender)

# 5 clusters
nmf_30.2 <- nmf.gene.expression(30,geneRpkm.2,zero.threshold = 0.5, log.scaled =TRUE)
nmf.30.2.rfcv.Dx <- rfcv(trainx=t(nmf_30.2$H[1:30,]),trainy=as.factor(pd.2$Dx),ntree=500,cv.fold=5,step=0.95)
nmf.30.2.rfcv.Race <- rfcv(trainx=t(nmf_30.2$H[1:30,]),trainy=as.factor(pd.2$Ethnicity),ntree=500,cv.fold=5,step=0.95)
nmf.boxplot(nmf_30.2$H,pd=pd.2,variable=pd.2$Gender,plot.layout = c(5,6))
nmf.25.2.rfcv.Gender <- rfcv(trainx=t(nmf_25.2$H[1:25,]),trainy=as.factor(pd.2$Gender),ntree=500,cv.fold=5,step=0.95)

par(mfrow=c(5,6))
for (i in 1:30){
  plot(pd.2$Age_of_Death,nmf_30.2$H[i,])
}

# 4 clusters
nmf_20.2 <- nmf.gene.expression(20,geneRpkm.2,zero.threshold = 0.01, log.scaled =TRUE)
nmf.20.2.rfcv.Dx <- rfcv(trainx=t(nmf_20.2$H[1:20,]),trainy=as.factor(pd.2$Dx),ntree=500,cv.fold=20,step=0.95)
nmf.20.2.rfcv.Race <- rfcv(trainx=t(nmf_20.2$H[1:20,]),trainy=as.factor(pd.2$Ethnicity),ntree=500,cv.fold=20,step=0.95)
nmf.boxplot(nmf_20.2$H,pd=pd.2,variable=pd.2$Dx,plot.layout = c(4,5))
nmf.20.2.rfcv.Gender <- rfcv(trainx=t(nmf_20.2$H[1:20,]),trainy=as.factor(pd.2$Gender),ntree=500,cv.fold=20,step=0.95)


######## experimenting with different loss types

nmf_15.1 <- nmf.gene.expression(15,geneRpkm.1,zero.threshold = 0.5, log.scaled = TRUE, loss.type = 'mse')
nmf.boxplot(nmf_15.1$H,pd=pd.1,variable=pd.1$Sex,plot.layout = c(3,5))
nmf.15.1.rfcv.Race <- rfcv(train=t(nmf_15.1$H[1:15,]),trainy=as.factor(pd.1$Race),ntree=500,cv.fold=5,step=0.95)
nmf.15.1.rfcv.Gender <- rfcv(train=t(nmf_15.1$H[1:15,]),trainy=as.factor(pd.1$Sex),ntree=500,cv.fold=5,step=0.95)


## 20 clusters?
nmf_20.1 <- nmf.gene.expression(20,geneRpkm.1,zero.threshold = 0.5, log.scaled = TRUE, loss.type = 'mse')
nmf.boxplot(nmf_20.1$H,pd=pd.1,variable=pd.1$Dx,plot.layout = c(4,5))
nmf.20.1.rfcv.Race <- rfcv(train=t(nmf_20.1$H[1:20,]),trainy=as.factor(pd.1$Race),ntree=500,cv.fold=5,step=0.95)
nmf.20.1.rfcv.Gender <- rfcv(train=t(nmf_20.1$H[1:20,]),trainy=as.factor(pd.1$Sex),ntree=500,cv.fold=5,step=0.95)
nmf.20.1.rfcv.Dx <- rfcv(train=t(nmf_20.1$H[1:20,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=5,step=0.95)


########### Decide on a k ########

# Choosing k is often the most important part of nnmf. By using missing value imputation, k can be chosen

geneRpkm.1.log <- log1p(geneRpkm.1)
par(mfrow=c(1,1))
set.seed(10)
plot(0, xlim = c(1,35), ylim = c(0,0.02), xlab = "Rank", ylab = "MSE")

k.values <- seq(10,20,5)
cols <- c('deepskyblue','orange','green','dark red')
i <- 1
err <- list()
for (i in 1:8) {
  
  sample.matrix <- sample(nrow(geneRpkm.1.log), size=1000)
  sample.geneRpkm.1.log <- geneRpkm.1.log[sample.matrix,]
  index <- sample(which(!is.na(geneRpkm.1.log[sample.matrix,])), 500);
  geneRpkm.temp <- sample.geneRpkm.1.log;
  geneRpkm.temp[index] <- NA;
  err[[i]] <- sapply(X = k.values,
                FUN = function(k, A) {
                  z <- nnmf(A, k, verbose = FALSE,n.threads=0);
                  mean((with(z, W%*%H)[index] - geneRpkm.1.log[sample.matrix,][index])^2)
                },
                A = geneRpkm.temp
  );
  #invisible(lines(err, col = col, type='b', lwd = 2, cex = 1));
  print(err[[i]])
  i = i + 1;
}


matrix.err <- matrix(unlist(err),ncol=3,byrow=TRUE)
plot(k.values,colMeans(matrix.err))

# investigate regularisation
beta.value <- 0.1
# first dataset, l1
nnmf_20.1.b1.0.1 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.1, zero.threshold = 0.1, 
                                        log.scaled = TRUE, 
                                        loss.type = "mse", b = c(0,0,beta.value))

nmf.20.1.b1.rfcv.Dx <- rfcv(train=t(nnmf_20.1.b1.0.1$H[1:20,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=10,step=0.8)
nmf.heatmap(nnmf_20.1.b1.0.1$H,variable=pd.1$Sex,step=2)

# first dataset, l2
nnmf_20.1.b2.0.1 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.1, zero.threshold = 0.1, 
                                        log.scaled = TRUE, 
                                        loss.type = "mse", b = c(beta.value,beta.value,0))

nmf.20.1.b2.rfcv.Dx <- rfcv(train=t(nnmf_20.1.b2.0.1$H[1:20,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=10,step=0.8)
nmf.heatmap(nnmf_20.1.b2.0.1$H,variable=pd.1$Dx,step=2)
print(nmf.20.1.b2.rfcv.Dx$error.cv)

# second dataset, l1
nnmf_20.2.b1.0.1 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.2, zero.threshold = 0.1, 
                                        log.scaled = TRUE, 
                                        loss.type = "mse", b = c(0,0,beta.value))
nmf.heatmap(nnmf_20.2.b1.0.1$H,variable=pd.2$Dx,step=2)

#second dataset, l2
nnmf_20.2.b2.0.1 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.2, zero.threshold = 0.1, 
                                        log.scaled = TRUE, 
                                        loss.type = "mse", b = c(beta.value,beta.value,0))





nmf.20.2.b1.rfcv.Dx <- rfcv(train=t(nnmf_20.2.b1.0.1$H[1:20,]),trainy=as.factor(pd.2$Dx),ntree=500,cv.fold=10,step=0.8)


nmf.heatmap <- function(H.matrix,variable,step=1){
  par(mar=c(1,1,1,1))
  H <- H.matrix/(rowMax(H.matrix)+1e-10)
  heatmap.2(H[,order(as.vector(variable))[seq(1,dim(H)[2],step)]],
            col=viridis::viridis(200),dendrogram='row',
            Rowv=TRUE,Colv=FALSE,trace='none',
            colsep=c(table(variable)[[1]]/step),
            sepwidth=c(1,0.05),
            labCol = FALSE,
            xlab = "Subjects",
            ylab = "Meta-Gene Expressions",
            main = "Meta-Gene Expression of Patients",
            density.info = "none",
            keysize = 1)
  par(mar=c(2,2,2,2))
}

nmf.heatmap(nnmf_20.2.b1.0.1$H,variable=pd.2$Dx,step=2)
nmf.boxplot(nnmf_20.1.b1.0.1$H,variable=pd.1$Race,pd=pd.1,plot.layout=c(4,5))








#### treating zeros as NA
nnmf_20.1.b1.0.1.5000 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.1, zero.threshold = 0.1, 
                                        log.scaled = TRUE, 
                                        loss.type = "mse", b = c(0,0,beta.value),
                                        na = FALSE,
                                        sample = TRUE,
                                        number.genes = 5000)
nmf.heatmap(nnmf_20.1.b1.0.1.500$H,variable=pd.1$Sex,step=2)
nmf.20.1.b1.na.rfcv.Dx <- rfcv(train=t(nnmf_20.1.b1.0.1.500$H[1:20,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=5,step=0.95)

nnmf_20.2.b1.0.1.5000 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.2, zero.threshold = 0.1, 
                                            log.scaled = TRUE, 
                                            loss.type = "mse", b = c(0,0,beta.value),
                                            na = FALSE,
                                            sample = TRUE,
                                            number.genes = 5000)
nmf.heatmap(nnmf_20.2.b1.0.1.5000$H,variable=pd.2$Gender,step=2)
nmf.20.2.b1.na.rfcv.Dx <- rfcv(train=t(nnmf_20.2.b1.0.1.5000$H[1:20,]),trainy=as.factor(pd.2$Dx),ntree=500,cv.fold=5,step=0.95)
print(nmf.20.2.b1.na.rfcv.Dx$error.cv)


### train leukemia

golub_train.nnmf <- nmf.gene.expression(n.clusters = 2,gene_expressions = golub_train.scaled, zero.threshold = 0.9,
                                        log.scaled = FALSE,
                                        loss.type = "mse", b = c(0,0,5),
                                        na = FALSE,
                                        sample = TRUE,
                                        number.genes = 1000)
nmf.heatmap(golub_train.nnmf$H,step=2,variable=as.character(golub_train$y))

nmf_10.1 <- nmf.gene.expression(n.clusters = 25, gene_expressions = geneRpkm.1, zero.threshold = 0.1, 
                                log.scaled = TRUE, 
                                loss.type = "mse",a = c(0.1,0.1,0), b = c(0.1,0.1,0),
                                na = FALSE,
                                sample = TRUE,
                                number.genes = 8000)
nmf.heatmap(nmf_10.1$H,step=1,variable=pd.1$Dx)
nmf.10.1.rfcv.Dx <- rfcv(train=t(nmf_10.1$H[1:25,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=10,step=0.95)
print(nmf.10.1.rfcv.Dx$error.cv)
nnmf.boxplot(nmf_10.1$H,variable = pd.1$Dx,plot.layout=c(5,5))







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

all_data.100.rfcv <- rfcv.gene.expression(folds=5,step=0.5,sample=TRUE,number.genes=2000,gene_expressions=geneRpkm.1,log.scaled=TRUE,zero.threshold=0.1,
                     variable = pd.1$Dx)


## finding responsible genes
nmf.responsible.genes <- function(H.matrix,W.matrix,variable,number.genes,gene.names,gene.expressions){
  fit=randomForest(y=as.factor(variable), x=t(H.matrix))
  importance.fit = importance(fit)
  print(importance.fit)
  responsible.meta.gene = which.max(importance.fit)
  top.gene_IDs <- (names(sort(W.matrix[,responsible.meta.gene],decreasing=TRUE)[1:number.genes]))
  top.gene.names <- gene.names[pmatch(top.gene_IDs,row.names(gene.expressions))]
  return(top.gene.names)
}

### heatmap for visualising gene expression datasets

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
pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//geneRpkm_1.pdf")
gene.expression.heatmap(geneRpkm.1,number.genes=500,log.scaled=T,row.scaled=T)
dev.off()

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//geneRpkm_2.pdf")
gene.expression.heatmap(geneRpkm.2,number.genes=500,log.scaled=T,row.scaled=T)
dev.off()

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//mice_image.pdf")
gene.expression.heatmap(mice_dat,number.genes=500,log.scaled=T,row.scaled=T,x.lab = "Mice")
dev.off()

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//golub_image.pdf")
gene.expression.heatmap(golub.positive,number.genes=500,log.scaled=T,row.scaled=T)
dev.off()

############ which genes are responsible ###############
consistent.genes <- list()
for (i in 1:5){
  
  nmf_25.1 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.1, zero.threshold = 0.1, 
                                  log.scaled = TRUE, 
                                  loss.type = "mse",a = c(0,0,0), b = c(0,0,0.1),
                                  na = FALSE,
                                  sample = TRUE,
                                  number.genes = 10000)
  
  nmf_25.2 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.2, zero.threshold = 0.1, 
                                  log.scaled = TRUE, 
                                  loss.type = "mse",a = c(0,0,0), b = c(0,0,0.1),
                                  na = FALSE,
                                  sample = TRUE,
                                  number.genes = 10000)
  
  #nmf.heatmap(nmf_25.1$H,step=1,variable=pd.1$Sex)
  #nmf.heatmap(nmf_25.2$H,step=1,variable=pd.2$Dx)
  
  
  
  genes.schizo.top.150 <- nmf.responsible.genes(nmf_25.1$H,nmf_25.1$W,variable=pd.1$Dx,number.genes=150,
                                                gene.names = geneMap$Symbol, gene.expressions = geneRpkm.1,
                                                control.age = T,age = pd.1$Age)
  
  
  genes.schizo.top.150.2 <- nmf.responsible.genes(nmf_25.2$H,nmf_25.2$W,variable=pd.2$Dx,number.genes=150,
                                                  gene.names = geneMap$Symbol, gene.expressions = geneRpkm.2,
                                                  control.age = T,age = pd.2$Age_of_Death)
  
  consistent.genes[[i]] <- genes.schizo.top.150.2[na.omit(pmatch(genes.schizo.top.150,genes.schizo.top.150.2))]
}
candidate.genes.schizo <- sort(table(unlist(consistent.genes)),decreasing=T)[1:7]

xtable(candidate.genes.schizo)

nmf_25.1 <- nmf.gene.expression(n.clusters = 20, gene_expressions = geneRpkm.1, zero.threshold = 0.1, 
                                log.scaled = TRUE, 
                                loss.type = "mse",a = c(0,0,0), b = c(0,0,0.1),
                                na = FALSE,
                                sample = TRUE,
                                number.genes = 10000)
nmf.heatmap(nmf_25.1$H,step=1,variable=pd.1$Sex)
nmf.heatmap(nmf_25.1$H,step=1,variable=pd.1$Race)
nmf.heatmap(nmf_25.1$H,step=1,variable=pd.1$Dx)

cor.dx.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Dx)))
cor.race.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Race)))
cor.age.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Age)))
cor.sex.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Sex)))
cor.df.2 <- as.data.frame(cbind(cor.dx,cor.race,cor.age,cor.sex))
names(cor.df) <- c("Schiz","Race","Age","Gender")
xtable(cor.df,caption="Correlation between metagenes from H and Class")



schizo.1.pca <- prcomp(t(pruned_expression_data))
schizo.pca.20 <- schizo.1.pca$x[,1:20]
cor.dx.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Dx)))
cor.race.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Race)))
cor.age.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Age)))
cor.sex.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Sex)))
cor.df.pca <- as.data.frame(cbind(cor.dx.pca,cor.race.pca,cor.age.pca,cor.sex.pca))
names(cor.df.pca) <- c("Schiz","Race","Age","Gender")
xtable(cor.df.pca,caption="Correlation between metagenes from PCA and Class")
genes.race <- nmf.responsible.genes(nmf_25.1$H,nmf_25.1$W,variable=pd.1$Race,number.genes=50,
                                    gene.names = geneMap$Symbol, gene.expressions = geneRpkm.1)

cor.dx.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Dx)))
cor.race.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Race)))
cor.age.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Age)))
cor.sex.2 <- cor(t(nmf_25.1$H),as.numeric(as.factor(pd.1$Sex)))
cor.df.2 <- as.data.frame(cbind(cor.dx,cor.race,cor.age,cor.sex))
names(cor.df) <- c("Schiz","Race","Age","Gender")
xtable(cor.df,caption="Correlation between metagenes from H and Class")



schizo.1.pca <- prcomp(t(pruned_expression_data))
schizo.pca.20 <- schizo.1.pca$x[,1:20]
cor.dx.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Dx)))
cor.race.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Race)))
cor.age.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Age)))
cor.sex.pca <- cor(schizo.pca.20,as.numeric(as.factor(pd.1$Sex)))
cor.df.pca <- as.data.frame(cbind(cor.dx.pca,cor.race.pca,cor.age.pca,cor.sex.pca))
names(cor.df.pca) <- c("Schiz","Race","Age","Gender")
xtable(cor.df.pca,caption="Correlation between metagenes from PCA and Class")
genes.race <- nmf.responsible.genes(nmf_25.1$H,nmf_25.1$W,variable=pd.1$Race,number.genes=50,
                                    gene.names = geneMap$Symbol, gene.expressions = geneRpkm.1)





pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//schizo_H_heatmap_1.pdf")
nmf.heatmap(nmf_25.1$H,step=3,variable=pd.1$Dx,variable.name="Disease")
dev.off()
pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//schizo_H_heatmap_2.pdf")
nmf.heatmap(nmf_25.2$H,step=3,variable=pd.2$Dx,variable.name="Disease")
dev.off()

nmf.20.1.rfcv.Dx <- rfcv(train=t(nmf_25.1$H[1:25,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=10,step=0.95)

nmf_20.1.b1.01.rfcv <- rfcv(train=t(nmf_25.1$H[1:20,]),trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=10,step=0.95)

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//random_forest_schizo_20.pdf")
par(mar=c(5,5,5,5))
par(mfrow=c(1,1))
plot(x=seq(20,1,-1),y=nmf_20.1.b1.01.rfcv$error.cv,ylab="error",
     xlab="dimensions used",main="Random Forrest Error against number of metagenes used",
     pch = 18)
dev.off()






################################## Additional portion on mice ##################################


load("mice_R_data.RData")


#### regressing out age

ge.1.corr.age <- cor(as.matrix(pd.1$Age), as.matrix(t(geneRpkm.1)))

ge.1.corr.age[is.na(ge.1.corr.age)] <- 0

ge.1.corr.sex <- cor(as.matrix(as.numeric(as.factor(pd.1$Sex))), as.matrix(t(geneRpkm.1)))

ge.1.corr.sex[is.na(ge.1.corr.sex)] <- 0

ge.1.corr.race <- cor(as.matrix(as.numeric(as.factor(pd.1$Race))), as.matrix(t(geneRpkm.1)))

ge.1.corr.race[is.na(ge.1.corr.race)] <- 0

ge.1.corr.dx <- cor(as.matrix(as.numeric(as.factor(pd.1$Dx))), as.matrix(t(geneRpkm.1)))

ge.1.corr.dx[is.na(ge.1.corr.dx)] <- 0

covariate.age.matrix <- ge.1.corr.age * t(ge.1.corr.age)

top.10.Dx.corr <- order(ge.1.corr.dx,decreasing=T)[1:10]
bottom.10.Dx.corr <- order(ge.1.corr.dx)[1:10]

rfcv.cor.dx <- rfcv(train=t(rbind(geneRpkm.1[.10.dx.corr,])),
                    trainy=as.factor(pd.1$Dx),ntree=500,cv.fold=10,step=0.95)


calculate_coph_corr <- function(k.values = seq(2,30,4),n.genes = 5000, zero.t = 0.1, n.runs = 5,gene_expressions){
  coph.cor <- list()
  k = 0
  for (i in k.values){
    k = k+1
    clusters.sum <- matrix(0,nrow = ncol(gene_expressions),ncol=i)
    for (j in 1:n.runs){
      
      nmf <- nmf.gene.expression(n.clusters = i, gene_expressions=gene_expressions,zero.threshold = zero.t,
                                   log.scaled = T, b=c(0,0,0),a=c(0,0,0),sample=T,number.genes=n.genes)
      H <- nmf$H
      clusters <- t(H)
      clusters[which(t(H)==rowMax(t(H)))] <- 1
      clusters[-which(t(H)==rowMax(t(H)))] <- 0
      clusters.sum <- clusters.sum + clusters
    }
    connectivity.matrix <- clusters.sum/n.runs
    d1 <- dist(connectivity.matrix)
    hc <- hclust(d1)
    d2 <- cophenetic(hc)
    coph.cor[[k]] <- c(i,cor(d1,d2))
    print(coph.cor[[k]])
  }
  
  coph_cors.vec <- unlist(coph.cor)
  coph_cors.mat <- cbind(coph_cors.vec[coph_cors.vec>1],coph_cors.vec[coph_cors.vec<=1])
  plot(coph_cors.mat)
  return(coph_cors.mat)
}


rfcv.nmf_12.1 <- rfcv(train=t(nmf_12.1$H),
                    trainy=as.factor(pd.1$Dx),ntree=1000,cv.fold=5,step=0.95)




######### benchmarking

benchmark(nmf.gene.expression(n.clusters = 5, gene_expressions = geneRpkm.1,log.scaled = T,zero.threshold = 0.1,
                              number.genes = 1000,a=c(0,0,0),b=c(0,0,0)),
          nmf.gene.expression(n.clusters = 5, gene_expressions = geneRpkm.1,log.scaled = T,zero.threshold = 0.1,
                              number.genes = 1000,a=c(0,0,0),b=c(0,0,0),use.base.nmf=T),
          replications = 0)




########### citations



########## initial dataset to test validity

data(golub)
golub.positive <- log1p(exp(golub))

golub.cl[golub.cl==0] <- "ALL"
golub.cl[golub.cl==1] <- "AML"

golub.nmf.2 <- nmf.gene.expression(n.clusters = 2,gene_expressions = golub.positive, zero.threshold = 0.1, log.scaled = F,
                    a = c(0,0,0), b = c(0,0,0.1))

golub.nmf.3 <- nmf.gene.expression(n.clusters = 3,gene_expressions = golub.positive, zero.threshold = 0.1, log.scaled = F,
                                   a = c(0,0,0), b = c(0,0,0.1))

golub.nmf.4 <- nmf.gene.expression(n.clusters = 4,gene_expressions = golub.positive, zero.threshold = 0.1, log.scaled = F,
                                   a = c(0,0,0), b = c(0,0,0.1))

golub.nmf.5 <- nmf.gene.expression(n.clusters = 5,gene_expressions = golub.positive, zero.threshold = 0.1, log.scaled = F,
                                   a = c(0,0,0), b = c(0,0,0.1))

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//golub_2.pdf")
nmf.heatmap(golub.nmf.2$H, variable = golub.cl, variable.name = "Tumor type")
dev.off()

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//golub_3.pdf")
nmf.heatmap(golub.nmf.3$H, variable = golub.cl, variable.name = "Tumor type")
dev.off()

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//golub_4.pdf")
nmf.heatmap(golub.nmf.4$H, variable = golub.cl, variable.name = "Tumor type")
dev.off()

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//golub_5.pdf")
nmf.heatmap(golub.nmf.5$H, variable = golub.cl, variable.name = "Tumor type")
dev.off()


### examining sparsity on the golub data set

nmf.heatmap(H.choice = F, W.choice = T, W.matrix = golub.nmf.4$W)


golub_coph_corrs <- calculate_coph_corr(n.genes = nrow(golub.positive), k.values = seq(2,5,1), zero.t = 0.1,
                                        n.runs = 10, gene_expressions = golub.positive)
old.par <- par(mar = c(0, 0, 0, 0))

par(mar = c(5,5,5,5))
plot.new()
pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//golub_coph_cor.pdf")
plot(golub_coph_corrs,xlab = "number of clusters", ylab = "cophenetic correlation coefficient",cex.lab = 2)

lines(golub_coph_corrs)
dev.off()
par(old.par)
########## mice dataset

mice.nmf.4 <- nmf.gene.expression(n.clusters = 6,gene_expressions = mice_dat, zero.threshold = 0.1, log.scaled = T,
                                   a = c(0,0,0), b = c(0,0,0),sample = T, number.genes = 5000)

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//aggression_nmf_H.pdf")
nmf.heatmap(mice.nmf.4$H, variable = aggression_labels, variable.name = "Mice Aggression")
dev.off()

aggression_genes <- nmf.responsible.genes(control.age = F, H.matrix = mice.nmf.4$H, W.matrix = mice.nmf.4$W, variable = aggression_labels,
                      number.genes = 20, gene.expressions = mice_dat,gene.names = row.names(mice_dat))



mice_coph_corrs <- calculate_coph_corr(n.genes = 5000, k.values = seq(2,8,1), zero.t = 0.1,
                                        n.runs = 10, gene_expressions = mice_dat)
old.par <- par(mar = c(0, 0, 0, 0))

par(mar = c(5,5,5,5))
plot.new()
pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//mice_coph_cor.pdf")
plot(mice_coph_corrs,xlab = "number of clusters", ylab = "cophenetic correlation coefficient",cex.lab = 2)

lines(mice_coph_corrs)
dev.off()
par(old.par)









############ Optional autoencoding?



##### Choosing k with value imputation

set.seed(567);

set.seed(10);
err <- list()
for (i in 1:5) {
  index2 <- sample(length(golub.positive), 2000);
  nsclc3 <- golub.positive;
  nsclc3[index2] <- NA;
  err[[i]] <- sapply(X = 2:8,
                FUN = function(k, A) {
                  z <- nnmf(A, k,verbose=F);
                  mean((with(z, W%*%H)[index2] - golub.positive[index2])^2)
                },
                A = nsclc3
  );

  
}
err.matrix.golub <- matrix(unlist(err),nrow=5,ncol=7,byrow = T)
par(mar=c(5,5,5,5))
pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//golub_imputations.pdf")
plot(seq(2,8,1),colMeans(err.matrix.golub), xlab = "clusters", ylab = "MSE", main = "Average MSE over 5 runs against the number of clusters used")

lines(seq(2,8,1),colMeans(err.matrix.golub))
dev.off()
gene_expressions = mice_dat
set.seed(10);
err.mice <- list()
pruned_expression_data <- gene_expressions


varying.row.index <- RowCV(pruned_expression_data)
pruned_expression_data <- pruned_expression_data[order(varying.row.index,decreasing=T)[1:5000],]
for (i in 1:5) {
  index3 <- sample(length(pruned_expression_data), 2000);
  nsclc4 <- pruned_expression_data;
  nsclc4[index3] <- NA;
  err.mice[[i]] <- sapply(X = 2:8,
                     FUN = function(k, A) {
                       z <- nnmf(A, k,verbose=T);
                       mean((with(z, W%*%H)[index3] - pruned_expression_data[index3])^2)
                     },
                     A = nsclc4
  );
}


err.matrix.mice <- matrix(unlist(err.mice),nrow=5,ncol=7,byrow = T)
par(mar=c(5,5,5,5))
pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//mice_imputations.pdf")
plot(seq(2,8,1),colMeans(err.matrix.mice), xlab = "clusters", ylab = "MSE", main = "Average MSE over 5 runs against the number of clusters used")

lines(seq(2,8,1),colMeans(err.matrix.mice))
dev.off()

gene_expressions = geneRpkm.1
set.seed(10);
err.schizo <- list()
zeros_geneRpkm <- gene_expressions==0

# how many zero entries per gene; look to reduce features
zeros_geneRpkm_per_gene <- rowSums(zeros_geneRpkm)

zero_gene_index <- which(zeros_geneRpkm_per_gene>0.1*dim(gene_expressions)[2])

#prune original gene data
print(paste(length(zero_gene_index), "dimensions removed"))
  pruned_expression_data <- gene_expressions[-zero_gene_index,]
  pruned_expression_data <- log1p(pruned_expression_data)

varying.row.index <- RowCV(pruned_expression_data)
pruned_expression_data <- pruned_expression_data[order(varying.row.index,decreasing=T)[1:10000],]
for (i in 1:3) {
  index4 <- sample(length(pruned_expression_data), 2000);
  nsclc5 <- pruned_expression_data;
  nsclc5[index4] <- NA;
  err.schizo[[i+6]] <- sapply(X = seq(28,32,4),
                          FUN = function(k, A) {
                            z <- nnmf(A, k,verbose=T,n.threads=0);
                            mean((with(z, W%*%H)[index4] - pruned_expression_data[index4])^2)
                          },
                          A = nsclc5
  );
}


err.matrix.schizo <- cbind(matrix(unlist(err.schizo[1:3]),nrow=3,ncol=6,byrow = T),matrix(unlist(err.schizo[7:9]),nrow=3,ncol=2,byrow = T))
par(mar=c(5,5,5,5))

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//schizo_imputations.pdf")
plot(seq(4,32,4),colMeans(err.matrix.schizo), xlab = "clusters", ylab = "MSE", main = "Average MSE over 3 runs against the number of clusters used")
lines(seq(4,32,4),colMeans(err.matrix.schizo))
dev.off()

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//schizo_coph_cors.pdf")
plot(coph_cors, xlab = "Number of Clusters", ylab = "cophenetic correlation coefficient")
lines(coph_cors)
dev.off()



##### classify using random forest

rfcv_nnmf_schizo.1 <- rfcv( trainx=t(nmf_25.1$H),trainy=as.factor(pd.1$Dx),ntree=1000,cv.fold=5,step=0.95)
rfcv_pca_schizo.1 <- rfcv( trainx=schizo.pca.20,trainy=as.factor(pd.1$Dx),ntree=1000,cv.fold=5,step=0.95)

rfcv_nnmf_race.1 <-  rfcv( trainx=t(nmf_25.1$H),trainy=as.factor(pd.1$Race),ntree=1000,cv.fold=5,step=0.95)
rfcv_pca_race.1 <- rfcv( trainx=schizo.pca.20,trainy=as.factor(pd.1$Race),ntree=1000,cv.fold=5,step=0.95)

rfcv_nnmf_gender.1 <- rfcv( trainx=t(nmf_25.1$H),trainy=as.factor(pd.1$Sex),ntree=1000,cv.fold=5,step=0.95)
rfcv_pca_gender.1 <- rfcv( trainx=schizo.pca.20,trainy=as.factor(pd.1$Sex),ntree=1000,cv.fold=5,step=0.95)

pdf("C://Users//aseem//OneDrive//dropbox//UCL_ML//Term_3//Project//Write-Up//Images//RF_plots.pdf")
plot.new()
par(mar=c(5,5,5,5))
plot(col="blue",type="l",ylim=c(0,0.4),rev(rfcv_nnmf_schizo.1$error.cv),xlab = "number of dimensions", 
     ylab = "5 fold cross validated error", main = "Random Forest Classifier on PCA and NMF generated features")
lines(col="blue",rev(rfcv_pca_schizo.1$error.cv),lty=2)
lines(col="red",rev(rfcv_nnmf_race.1$error.cv))
lines(col="red",rev(rfcv_pca_race.1$error.cv),lty=2)
lines(col="dark green",rev(rfcv_nnmf_gender.1$error.cv))
lines(col=" dark green",rev(rfcv_pca_gender.1$error.cv),lty=2)
legend(col=c("blue","red","dark green"),legend=c("PCA","NMF","PCA"))
dev.off()





