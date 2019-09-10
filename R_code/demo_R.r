# Copyright (C) 2019 The University of Edinburgh 
# Author Stuart Aitken MRC IGMM stuart.aitken@igmm.ed.ac.uk
# All Rights Reserved.
# Funded by the Medical Research Council
# https://www.ed.ac.uk/mrc-human-genetics-unit


library(e1071);
library(Hmisc);

source('functions.r');

quartz(width=12,height=4)
par(mfrow=c(1,2));

# example 1. 5 unlabelled data points
# 5 data points (4 discrete attributes) with two distinct patterns
obs <- rbind(c( "1", "1","-1","-1"),     #  pattern A
             c("-1","-1", "1", "1"),     #  pattern B
             c( "1", "1","-1","-1"),     #  pattern A
             c("-1","-1", "1", "1"),     #  pattern B
             c("-1","-1", "1", "1"));    #  pattern B
colnames(obs) <- c('d1','d2','d3','d4');
obs <- data.frame(obs);

# prior to calling EM Naive Bayes, randomly initialise probability tables following the naive Bayes object convention
# q0 specifies the initial class probabilities
# q_y0 is the list of attribute values / class probability tables
set.seed(2017);
q0   <- c(x<-runif(1),1-x);
tables2 <- vector('list',length=4);
for(i in 1:4) {
    tables2i <- rbind(c(x<-runif(1),1-x),
                       c(x<-runif(1),1-x));
    tables2i <- data.frame(tables2i);
    colnames(tables2i) <- c('-1','1');
    rownames(tables2i) <- c('y1','y2');
    tables2[[i]] <- tables2i;
}
tables2;

# call EM Naive Bayes specifying 2 class labels, plot log likelihood on completion of EM
resultEM_k_eq_2 <- EM_NB_Fn(obs,q0,tables2,maxiter=10,classLabels=c('y1','y2'),PLOT=TRUE);
resultEM_k_eq_2;

# predict labels for data using standard predict() method
pred_class <- predict(as_NBobj(resultEM_k_eq_2,colnames(obs)),obs);
pred_class;



# example 2. iris data
data(iris);

# discretise iris data into 3 bins
irisDiscr <- iris;
irisDiscr$Sepal.Length <- cut2(irisDiscr$Sepal.Length,g=3);
irisDiscr$Sepal.Width  <- cut2(irisDiscr$Sepal.Width,g=3);
irisDiscr$Petal.Length <- cut2(irisDiscr$Petal.Length,g=3);
irisDiscr$Petal.Width  <- cut2(irisDiscr$Petal.Width,g=3);
irisDiscr$Species <- factor(irisDiscr$Species);

# prior to calling EM Naive Bayes, randomly initialise probability tables for 3 classes 
# q01, q02, (1-q01-q02)  specify the initial class probability
# tables3 is the list of attribute values / class probability tables [for simplicity here, assume 3 bins in all cases]
set.seed(2109);

q01 <- runif(1,min=0,max=0.4);
q02 <- runif(1,min=0,max=0.4);
tables3 <- vector('list',length=4);
for(i in 1:4) {
    tables3[[i]] <- (array(dim=c(3,3),0));
    rownames(tables3[[i]]) <- c('y1','y2','y3');
    colnames(tables3[[i]]) <- levels(irisDiscr[,i]);
    for(j in 1:3) {
        q01ij <- runif(1,min=0,max=0.4);
        q02ij <- runif(1,min=0,max=0.4);
        tables3[[i]][j,] <- c(q01ij,q02ij,1-(q01ij+q02ij));
    }
}
tables3;

# call EM Naive Bayes specifying 3 class labels, plot log likelihood on iteration 20 and on completion of EM
resultEM_k_eq_3 <- EM_NB_Fn(irisDiscr[,1:4],c(q01,q02,(1-q01-q02)),tables3,
                            maxiter=100,classLabels=c('y1','y2','y3'),FACTORS=TRUE,PLOT=TRUE);
lll <- length(resultEM_k_eq_3$ll);

cat(paste('no. iterations:',(lll-1),'last delta LL:',resultEM_k_eq_3$ll[lll]-resultEM_k_eq_3$ll[(lll-1)],'\n'));


# predict labels for data using standard predict() method
pred_class_iris <- predict(as_NBobj(resultEM_k_eq_3,colnames(irisDiscr)[1:4]),irisDiscr[,1:4]);
pred_class_iris;

# compare with true labels
table(pred=pred_class_iris,true=irisDiscr$Species);

#pred setosa versicolor virginica
#  y3     50          0         0
#  y2      0         49        10
#  y1      0          1        40

par(mfrow=c(1,3)); # compare with PCA projection
x <- prcomp(iris[,1:4]);
cols <- array(dim=length(pred_class_iris),'red')
cols[pred_class_iris=='y2']<-'green'
cols[pred_class_iris=='y3']<-'blue'
plotPCA(x,choices=c(1,2),cols);
plotPCA(x,choices=c(1,3),cols);
plotPCA(x,choices=c(1,4),cols);



# example 3. iris data looking for 4 classes

# prior to calling EM Naive Bayes, randomly initialise probability tables for 4 classes 
set.seed(2099);
q01 <- runif(1,min=0,max=0.3);
q02 <- runif(1,min=0,max=0.3);
q03 <- runif(1,min=0,max=0.3);
tables4 <- vector('list',length=4);
for(i in 1:4) {
    tables4[[i]] <- (array(dim=c(4,3),0));
    rownames(tables4[[i]]) <- c('y1','y2','y3','y4');
    colnames(tables4[[i]]) <- levels(irisDiscr[,i]);
    for(j in 1:4) {
        q01ij <- runif(1,min=0,max=0.3);
        q02ij <- runif(1,min=0,max=0.3);
        tables4[[i]][j,] <- c(q01ij,q02ij,1-(q01ij+q02ij));
    }
}
tables4;

par(mfrow=c(1,1)); 
resultEM_k_eq_4 <- EM_NB_Fn(irisDiscr[,1:4],c(q01,q02,q03,(1-q01-q02-q03)),tables4,
                            maxiter=500,classLabels=c('y1','y2','y3','y4'),FACTORS=TRUE,PLOT=TRUE);
lll <- length(resultEM_k_eq_4$ll);

cat(paste('no. iterations:',(lll-1),'last delta LL:',resultEM_k_eq_4$ll[lll]-resultEM_k_eq_4$ll[(lll-1)],'\n'));


# predict labels for data using standard predict() method
pred_class_iris_4 <- predict(as_NBobj(resultEM_k_eq_4,colnames(irisDiscr)[1:4]),irisDiscr[,1:4]);
pred_class_iris_4;

# compare with true labels
table(pred=pred_class_iris_4,true=irisDiscr$Species);

#pred setosa versicolor virginica
#  y2     50          0         0
#  y1      0         49         9
#  y3      0          0        23
#  y4      0          1        18

