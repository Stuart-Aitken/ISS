# Copyright (C) 2019 The University of Edinburgh 
# Author Stuart Aitken MRC IGMM stuart.aitken@igmm.ed.ac.uk
# All Rights Reserved.
# Funded by the Medical Research Council
# https://www.ed.ac.uk/mrc-human-genetics-unit


DBL.MAX <- 1.79769e+308;

productOverAttributes <- function(l,y,x) {
    r <- 1;
    for(j in 1:length(l)) { # j in 1..d attributes
        if(!is.na(x[j])) {
            r <- r*l[[j]][y,as.character(x[j])];
        }
    }
    return(r);
}
productOverAttributesForFactors <- function(l,y,x) {
    r <- 1;
    for(j in 1:length(l)) { # j in 1..d attributes
        if(!is.na(x[j])) {
            r <- r*l[[j]][y,unlist(x[j])];
        }
    }
    return(r);
}

checkSum1Fn <- function(q_y0) {
    x <- 0;
    for(i in 1:length(q_y0)) {
        for(j in 1:length(q_y0[[i]][,1])) {
            if(abs(1-sum(q_y0[[i]][j,]))>5e-16) { x<-x+1; cat(paste(i,j,abs(1-sum(q_y0[[i]][j,])),'\n')); print(q_y0[[i]][j,]); }
        }
    }
    if(x>0) { print(paste('no 1 errors:',x)); }
}


log.plus <- function(x,y) {
    if(x>y) x+log(1+exp(y-x))
    else    y+log(1+exp(x-y))
}

loglikeFn <- function(obs,q0,q_y0,classes=NULL,FACTORS=FALSE,llForPZero=log(1e-6),REPORT=FALSE) {
    if(!is.null(classes) && sum((classes))==prod(dim(classes))) {
        foundNA <- FALSE;
        llFn <- function(obsi) {
            sy <- -DBL.MAX;                      # log of sum likelihood for all [applicable] classes
            for(y in 1:length(q_y0[[1]][,1])) {  # sum for all classes
                if(classes[i,y]) {               # test if y allowable for xi
                    piy <- log(q0[y]);
                    for(j in 1:length(q_y0)) {   # product for all attributes [sum of log]
                        if(!is.na(obsi[j])) {
                            if(FACTORS) {
                                piy <- piy + max(llForPZero,log(q_y0[[j]][y,unlist(obsi[j])]));  # when p=0 log(0)=-Inf
                            } else {
                                piy <- piy + max(llForPZero,log(q_y0[[j]][y,as.character(obsi[j])]));
                            }
                        } else {
                            piy <- piy + log(1/length(q_y0[[j]][y,])); #   if NA use p = 1/no. values
                            foundNA <- TRUE;
                        }
                    } # piy is the sum log(q0[y]) + log(q_y0) for all attributes
                    sy <- log.plus(sy,piy); # add log values sy and piy
                }
            }
            return(sy);
        }
        ll <- array(apply(obs,1,llFn));
    } else {
        if(is.null(classes)) {
            classes <- array(dim=c(length(obs[,1]),length(q_y0[[1]][,1])),TRUE);
        }
        ll <- array(dim=length(obs[,1]),NA);
        foundNA <- FALSE;
        for(i in 1:length(obs[,1])) {            # sum for all obs
            sy <- -DBL.MAX;                      # log of sum likelihood for all [applicable] classes
            for(y in 1:length(q_y0[[1]][,1])) {  # sum for all classes
                if(classes[i,y]) {               # test if y allowable for xi
                    piy <- log(q0[y]);
                    for(j in 1:length(q_y0)) {   # product for all attributes [sum of log]
                        if(!is.na(obs[i,j])) {
                            if(FACTORS) {
                                piy <- piy + max(llForPZero,log(q_y0[[j]][y,unlist(obs[i,j])]));  # when p=0 log(0)=-Inf
                            } else {
                                piy <- piy + max(llForPZero,log(q_y0[[j]][y,as.character(obs[i,j])]));
                            }
                        } else {
                            piy <- piy + log(1/length(q_y0[[j]][y,])); #   if NA use p = 1/no. values
                            foundNA <- TRUE;
                        }
                    } # piy is the sum log(q0[y]) + log(q_y0) for all attributes
                    sy <- log.plus(sy,piy); # add log values sy and piy
                }
            }
            ll[i] <- sy;
        }
    }
    if(REPORT && foundNA) { cat('found NAs in loglikelihood calculation\n'); }
    if(REPORT && sum(ll==-Inf)>0) {
        cat('found -Inf in loglikelihood calculation\n'); 
    }
    return(sum(ll[ll!=-Inf]));
}


nonDecreasingFn <- function(a) {
    r <- TRUE;
    if(length(a)>1) {
        for(i in 2:length(a)) {
            if(a[i-1]>a[i]) { cat(paste(a[i-1]-a[i]),';',sep=''); r <- FALSE; }
        }
    }
    return(r);
}

EM_NB_Fn <- function(obs,q0,q_y0,maxiter=10,minDeltaLL=0.000001,classLabels=c('y1','y2'),FACTORS=FALSE,PLOT=FALSE,classes=NULL) {
    if(!is.null(classes) && sum((classes))==prod(dim(classes))) { cat('speeding up...\n'); }
    if(!is.null(classes)) {
        classLabels <- colnames(classes);
        samplesClassCount <- apply(classes,1,sum);
        if(sum(samplesClassCount==0)>0) { cat(paste('** error: ',sum(samplesClassCount==0),'samples with no possible class label\n')); return(NULL); }
    } else {
        classes <- array(dim=c(length(obs[,1]),length(classLabels)),TRUE);
        colnames(classes) <- classLabels;
    }
    ll <- loglikeFn(obs,q0,q_y0,classes=classes,FACTORS=FACTORS,REPORT=F);
    #cat(paste('logLikelihood',signif(ll,6),'\n'));
    noObs <- length(obs[,1]);
    noAttr <- length(obs[1,]);
    noClasses <- length(classLabels);
    noClassesForObs <- apply(classes,1,sum);
    noObsForClasses <- apply(classes,2,sum);
    expectedCountPerClass <- array(dim=noClasses);
    rownames(expectedCountPerClass) <- classLabels;
    for(i in 1:noClasses) {
        expectedCountPerClass[i] <- sum(1/noClassesForObs[classes[,i]]);
    }
    for(iter in 1:maxiter) {
        #cat(paste(iter,(proc.time()-t0)[1],';'));  t0 <- proc.time();
        delta <- array(dim=c(noClasses,noObs),NA);  # the 'responsibility' array: classes (rows) samples (cols)
        colnames(delta) <- rownames(obs);
        rownames(delta) <- classLabels;    # get delta y i [without dividing by sum over i]
        if(!is.null(classes) && sum((classes))==prod(dim(classes))) {
            if(FACTORS) {
                for(y in 1:noClasses) {
                    delta[y,] <- apply(obs,1,function(obsi) { q0[y]*productOverAttributesForFactors(q_y0,y,obsi); });
                }
            } else {
                for(y in 1:noClasses) {
                    delta[y,] <- apply(obs,1,function(obsi) { q0[y]*productOverAttributes(q_y0,y,obsi); });
                }
            }
        } else {
            for(i in 1:noObs) {                # i in 1..n samples
                for(y in 1:noClasses) {        # j in 1..k classes
                    if(FACTORS) {
                        if(classes[i,y]) {
                            delta[y,i] <- q0[y]*productOverAttributesForFactors(q_y0,y,obs[i,]);
                        } else {
                            delta[y,i] <- NA;
                        }
                    } else {
                        if(classes[i,y]) {
                            delta[y,i] <- q0[y]*productOverAttributes(q_y0,y,obs[i,]);
                        } else {
                            delta[y,i] <- NA;
                        }
                    }
                }
            }
        }
        sumYs <- apply(delta,2,sum,na.rm=T);       # get sum over i and divide
        
        for(i in 1:noObs) {
            if(sumYs[i]!=0) {
                delta[,i] <- delta[,i]/sumYs[i];
            } else {
                delta[classes[i,],i] <- min(sumYs[sumYs>0])/noClassesForObs[i];  # set values for rows in delta that can be assigned values
            }
        }

        # update q0 to qt [sum deltas in row / no obs]
        
        qt <- apply(delta,1,sum,na.rm=T)/noObs; # i in 1..n samples  *** note possible differing no obs / class
        
        if(sum(qt)!=1) { qt <- qt/sum(qt); }
        # update q_y0 to q_yt 
        q_yt <- vector('list',length=noAttr);
        for(i in 1:noAttr) {
            q_yti <- NULL;
            iValueSet <- sort(unique(obs[,i]));
            noValues <- length(iValueSet);
            for(y in 1:noClasses) {
                q_yti_vals <- NULL;
                for(valIndx in 1:noValues) {
                    q_yti_vals <- c(q_yti_vals,    # note selection for classes[,y] allowable classes
                                    sum(delta[y,obs[,i]==iValueSet[valIndx] & classes[,y]],na.rm=T)/sum(delta[y,classes[,y]],na.rm=T));
                }
                if(sum(q_yti_vals)==0) {
                    q_yti_vals <- array(dim=noValues,1/noValues);
                } else if(sum(q_yti_vals)!=1) {
                    q_yti_vals <- q_yti_vals/sum(q_yti_vals);
                }
                q_yti <- rbind(q_yti,q_yti_vals);
            }
            q_yti <- data.frame(q_yti,row.names=classLabels);
            colnames(q_yti) <- as.character(iValueSet);
            q_yt[[i]] <- q_yti;
        }
        #print(qt);
        #print(q_yt);

        checkSum1Fn(q_yt);
        
        q0 <- qt;
        q_y0 <- q_yt;
        lli <- loglikeFn(obs,q0,q_y0,classes=classes,FACTORS=FACTORS);
        ll <- c(ll,lli);

        if(iter>(maxiter/3) && ((!is.null(minDeltaLL) && abs(ll[iter]-ll[(iter-1)])<=minDeltaLL) || ll[iter]==ll[iter-1])) {
            if(PLOT) {
                stindx <- seq(from=1,to=(iter+1))[ll>(lli-5)][1]-1; 
                plot(ll[1:length(ll)],type='l',main=paste('Iteration',iter),xlab='Iteration',ylab='Log likelihood');
            }
            if(!nonDecreasingFn(ll)) { cat('\nWARNING: log L decreases\n'); }
            return(list(qt=qt,q_yt=q_yt,ll=ll)); }
    }

    
    if(!nonDecreasingFn(ll)) { cat('\nWARNING: log L decreases\n'); }
    return(list(qt=qt,q_yt=q_yt,ll=ll));
}


# create a "naiveBayes" object and set its probability tables to those of EM_NB_Fn() output
as_NBobj <- function(emo,attrnames) {
    noattrs <- unlist(lapply(emo$q_yt,function(x) { dim(x)[2];}));
    noclasses <- unique(unlist(lapply(emo$q_yt,function(x) { dim(x)[1];})));
    max_noattrs <- max(noattrs,noclasses);
    obs <- data.frame(array(dim=c(max_noattrs,length(noattrs))));
    for(i in 1:length(noattrs)) {
        obs[,i] <- as.factor(array(dim=max_noattrs,colnames(emo$q_yt[[i]])));
    }
    colnames(obs) <- attrnames;
    y <- as.factor(array(dim=max_noattrs,rownames(emo$q_yt[[1]])));
    nbobj <- naiveBayes(y ~ ., data=data.frame(obs,y));
    
    for(i in 1:length(nbobj$tables)) {
        for(j in 1:length(nbobj$tables[[i]][1,])) {
            for(k in 1:length(nbobj$tables[[i]][,1])) {
                nbobj$tables[[i]][k,j] <- emo$q_yt[[i]][k,j];
            }
        }
    }
    nbobj$apriori  <- (emo$qt * max_noattrs);
    return(nbobj);
}


plotPCA <- function(pcres,choices=1L:2L,cols=NULL) {
    if(is.null(cols)) { cols <- array(dim=dim(pcres$x)[1],'red'); }
    scores <- pcres$x
    lam <- pcres$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    scale <- 1;
    if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
    if(scale != 0) lam <- lam^scale else lam <- 1
    pts<-t(t(scores[, choices]) / lam);
    
    plot(pts,col=cols, pch=19,main='',lwd=2);
}
