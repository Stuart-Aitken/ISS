# Copyright (C) 2019 The University of Edinburgh 
# Author Stuart Aitken MRC IGMM stuart.aitken@igmm.ed.ac.uk
# All Rights Reserved.
# Funded by the Medical Research Council
# https://www.ed.ac.uk/mrc-human-genetics-unit


# EM Naive Bayes in Java
# The following command will apply the Java EM algorithm 1000 times for each of k=2..10 to
# the 150 samples of discretised data in iris.csv starting each run of 1000 with the random seed 101
# [the best solution for each value of k is retained]:
# ./run.sh 101 iris.csv 

# For each value of k, text files with the model probabilities
# and with the trace of log likelihood values are saved,
# for k=2:
# EMNB_101_2_iris_model.txt
# EMNB_101_2_iris_logls.txt

# The results can be read with the R code below [a tar file of precomputed results is also provided: Results.tar.gz].
# Note that the number of values per attribute is fixed at 3 for 4 data attributes in NBO.java:
#    public static int[] novalues = {3,3,3,3}; 
# and that the values in the discretised data file must be {-1,1,NA} or {-1,0,1,NA}
# (these are fixed in DO.java to map to array index values 0,1 or 0,1,2; -1 is used for NA).
# The number of values per attribute is also fixed in the R code below,
# but the code can be generalised in the obvious way.


library(e1071);
library(Hmisc);


# create a "naiveBayes" object and set its probability tables to those of the Java output
fileToNBayesObject <- function(noY,obs,fni) {
    java_model <- read.table(fni,header=F,sep=','); 
    noValues <- c(3,3,3,3); 
    startColInFile <- cumsum(c(2,3,3,3)); 
    classLabels <- paste('y',seq(1:noY),sep='');
    y <- as.factor(array(dim=dim(obs)[1],classLabels));
    obsy <- data.frame(obs,y,stringsAsFactors=T);
    
    nb <- naiveBayes(y ~ ., data = obsy);
    nb$apriori <- java_model[,1]*length(obs[,1]);
    for(i in 1:length(nb$tables)) {
        if(noValues[i]==2 && length(nb$tables[[i]][1,])==2) {
           nb$tables[[i]][,1] <- java_model[,startColInFile[i]];
           nb$tables[[i]][,2] <- java_model[,startColInFile[i]+1];
        } else if(noValues[i]==3 && length(nb$tables[[i]][1,])==3) {
           nb$tables[[i]][,1] <- java_model[,startColInFile[i]];
           nb$tables[[i]][,2] <- java_model[,startColInFile[i]+1];
           nb$tables[[i]][,3] <- java_model[,startColInFile[i]+2];           
        } else { print('error'); }
    }
    names(nb$tables)<- colnames(obs);
    return(nb);
}

# the seed used in the Java run was 101, this forms part of the file name
seed <- 101;

# read each set of results for each value of k
iris <- read.table(file='iris.csv',sep=',',header=T,stringsAsFactors=T);
iris[,1] <- as.factor(iris[,1]);
iris[,2] <- as.factor(iris[,2]);
iris[,3] <- as.factor(iris[,3]);
iris[,4] <- as.factor(iris[,4]);
resultListModels <- vector('list',length=10);
resultListLogL <- vector('list',length=10);

for(i in 2:10) {
    fni <- paste('EMNB_',seed,'_',i,'_iris_model.txt',sep='')
    fnlli <- paste('EMNB_',seed,'_',i,'_iris_logls.txt',sep='')
    if(file.exists(fni)) {
        resultListModels[[i]] <- fileToNBayesObject(i,iris[,1:4],fni);
        lli <- read.table(fnlli,header=F,sep=',')[,1]; 
        resultListLogL[[i]] <-  rev(lli[lli!=0])[1];
    } else {
        cat(paste('warning no results for',i,'\n'));
    }
}


# the final log likelihood values
logLs <- unlist(resultListLogL);
# 

# calculate AIC accounting for free parameters: 8 attributes have 2 free parameters; 2 have 1 free parameter giving
# (4*2) plus (4*2+1) * no. classes beyond 1; the additional +1 is for the probabilities of the k classes
aics <- -2*logLs + 2*((4*2) + (4*2+1)*(1:9));
# 

indxMin <- seq(from=2,to=10)[aics==min(aics)];
indxMin;
# k=3 is optimal


plot(2:10,aics,ylab='AIC',xlab='No. classes',pch=19);
points(indxMin,aics[indxMin-1],col='red',pch=19,cex=1.1);

resultListModels[[indxMin]];

pred_class_iris <- predict(resultListModels[[indxMin]],iris[,1:4])


# compare with true labels
table(pred=pred_class_iris,true=iris$Species);

#    true
#pred  1  2  3
#  y3 50  0  0
#  y1  0 49 10
#  y2  0  1 40

