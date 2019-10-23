
DRRCI_counts <- function(x1, n1, x2, n2, x0, n0, x0b=NULL, n0b=NULL, options=NULL) {

  op  <- check.options(options)
  tmp <- checkCounts(x1, n1, x2, n2, x0, n0, x0b, n0b, op)

  ret <- DRRCI_main(x1=tmp$x1,x2=tmp$x2,x3=tmp$x3,x4=tmp$x4,
                    n1=tmp$n1,n2=tmp$n2,n3=tmp$n3,n4=tmp$n4,
                    tr=NA,Y=NA,Z=NA,
                    alpha=op$alpha, altParam=op$altParam, deltaMethod=op$deltaMethod,
                    refPop=op$refPop, LRT=op$LRT, estRR=0, op=op) 

  ret

} # END: DRRCI_counts

DRRCI_unadj <- function(trtGroup, outcome, strata=NULL, options=NULL) {

  op       <- check.options(options)
  checkData(trtGroup, outcome, S=strata)
  tmp      <- setupData(trtGroup, outcome, op, S=strata)

  ret <- DRRCI_main(x1=tmp$x1,x2=tmp$x2,x3=tmp$x3,x4=tmp$x4,
                    n1=tmp$n1,n2=tmp$n2,n3=tmp$n3,n4=tmp$n4,
                    tr=NA,Y=NA,Z=NA,
                    alpha=op$alpha, altParam=op$altParam, deltaMethod=op$deltaMethod,
                    refPop=op$refPop, LRT=op$LRT, estRR=0, op=op) 

  ret

} # END: DRRCI_unadj


DRRCI_adj <- function(trtGroup, outcome, covars, options=NULL) {

  if (!length(covars)) stop("ERROR: covars must be specified with DRRCI_adj")
  op       <- check.options(options, method="adj")
  checkData(trtGroup, outcome, X=covars)
  tmp      <- setupData(trtGroup, outcome, op, X=covars)
  trtGroup <- tmp$trt
  outcome  <- tmp$Y
  covars   <- tmp[["X", exact=TRUE]]
  rm(tmp)
  gc()

  ret <- DRRCI_main(x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,
               tr=trtGroup,Y=outcome,Z=covars,
               alpha=op$alpha, altParam=0, deltaMethod=0,
               refPop="NULL", LRT=op$LRT, estRR=op$estRR, op=op) 

  ret

} # END: DRRCI_adj

setProbs <- function(vec, MINPROB=1e-6, MAXPROB=1-1e-6, stopOnError=1) {

  x <- vec

  # Check for errors
  if (stopOnError) {
    if (any(vec < 0)) stop("ERROR with probabilities")
    if (any(vec > 1)) stop("ERROR with probabilities")
  }

  tmp                      <- x < MINPROB
  tmp2                     <- is.na(tmp)
  if (any(tmp2)) tmp[tmp2] <- FALSE
  if (any(tmp)) x[tmp]     <- MINPROB
  tmp                      <- x > MAXPROB
  tmp2                     <- is.na(tmp)
  if (any(tmp2)) tmp[tmp2] <- FALSE
  if (any(tmp)) x[tmp]     <- MAXPROB

  x

} # END: setProbs

setupX <- function(X, nrows) {

  cx <- NULL
  if (length(X)) cx <- colnames(X)  

  if (!length(X)) {
    ret <- matrix(data=1, nrow=nrows, ncol=1)
    colnames(ret) <- "Intercept"
    return(ret)
  }

  m    <- ncol(X)
  keep <- rep(TRUE, m)
  if (!length(cx)) cx <- paste("Covar", 1:m, sep="")
  for (i in 1:m) {
    tmp <- var(X[, i])
    if (!is.finite(tmp) || (tmp == 0)) keep[i] <- FALSE
  }
  if (!all(keep)) {
    X  <- X[, keep, drop=FALSE]
    cx <- cx[keep]
  }
  if (!ncol(X)) {
    X  <- matrix(data=1, nrow=nrows, ncol=1)
    cx <- "Intercept"
  } else {
    X <- cbind(1, X)
    cx <- c("Intercept", cx)
  }
  colnames(X) <- cx

  X

} # END: setupX

getRowsToUse <- function(trt, Y, X, S) {

  ret <- is.finite(Y) & is.finite(trt) 
  if (length(X)) ret <- ret & (!rowSums(!is.finite(X)))
  if (length(S)) ret <- ret & is.finite(S)

  ret

} # END: getRowsToUse

setupData <- function(trt, Y, op, X=NULL, S=NULL) {

  if (!is.vector(trt)) trt <- makeVector(trt)
  if (!is.vector(Y)) Y <- makeVector(Y)
  if (length(X) && !is.matrix(X)) X <- as.matrix(X)
  if (length(S)) S <- makeVector(S)

  tmp <- getRowsToUse(trt, Y, X, S)
  if (!all(tmp)) {
    Y   <- Y[tmp]
    trt <- trt[tmp]
    if (length(X)) X <- X[tmp, , drop=FALSE]
    if (length(S)) S <- S[tmp]
  }
  if (!length(Y)) stop("ERROR: all observations have been removed")

  tmp <- !(Y %in% 0:1)
  if (any(tmp)) stop("ERROR: outcome must be coded 0 or 1")
  utrt <- sort(unique(trt))
  ntrt <- length(utrt)
  if (length(X) && (ntrt != 3)) stop("ERROR: there can only be 3 treatment groups for this method")
  if ((ntrt != 3) && (ntrt != 4)) stop("ERROR: there must be 3 or 4 treatment groups")
  if ((ntrt == 3) && !all(utrt == 0:2)) stop("ERROR: trtGroup must be coded 0-2")
  if ((ntrt == 4) && !all(utrt == 1:4)) stop("ERROR: trtGroup must be coded 1-4")
  #if (ntrt == 4) trt <- trt - 1
    
  if (length(X)) {
    X   <- setupX(X, length(Y))
    ret <- list(Y=Y, trt=trt, X=X)
  } else {
    ret <- getCounts(trt, Y, S, ntrt, op)
  }

  ret

} # END: setupData

checkCounts <- function(x1, n1, x2, n2, x3, n3, x4, n4, op) {

  #  x3, x4 are control pops

  checkVector(x1, "x1", len=0)
  len <- length(x1)
  checkVector(n1, "n1", len=len)
  checkVector(x2, "x2", len=len)
  checkVector(n2, "n2", len=len)
  checkVector(x3, "x3", len=len)
  checkVector(n3, "n3", len=len)
  if (op$contCorrect) {
    x1[x1==0] <- 0.5
    n1[n1==0] <- 0.5
    x2[x2==0] <- 0.5
    n2[n2==0] <- 0.5
    x3[x3==0] <- 0.5
    n3[n3==0] <- 0.5
  }
  if (any(n1 < x1)) stop("ERROR with x1 and/or n1")
  if (any(n2 < x2)) stop("ERROR with x2 and/or n2")
  if (any(n3 < x3)) stop("ERROR with x3 and/or n3")
  if (length(x4) || length(n4)) {
    checkVector(x4, "x4", len=len)
    checkVector(n4, "n4", len=len)
    if (op$contCorrect) {
      x4[x4==0] <- 0.5
      n4[n4==0] <- 0.5
    }
    if (any(n4 < x4)) stop("ERROR with x4 and/or n4")
    names(x4) <- NULL
    names(n4) <- NULL 
  } else {
    x4 <- NA
    n4 <- NA
  }
  
  names(x1) <- NULL
  names(x2) <- NULL
  names(x3) <- NULL
  names(n1) <- NULL
  names(n2) <- NULL
  names(n3) <- NULL

  list(x1=x1, n1=n1, x2=x2, n2=n2, x3=x3, n3=n3, x4=x4, n4=n4)

} # END: checkCounts

getCountVecs <- function(utrt, trt, Y1, S, uS, nS, op) {

  x    <- rep(NA, nS)
  n    <- x

  tmpT <- trt == utrt
  if (any(tmpT)) {
    for (j in 1:nS) {
      tmpS <- S == uS[j]
      x[j] <- sum(Y1 & tmpT & tmpS)
      n[j] <- sum(tmpT & tmpS)
    } 
    if (op$contCorrect) {
      tmp <- n == 0
      if (any(tmp)) n[tmp] <- 0.5
      tmp <- x == 0
      if (any(tmp)) x[tmp] <- 0.5
    }
  } 

  list(x=x, n=n)

} # END: getCountVec

getCounts <- function(trt, Y, S, ntrt, op) {

  # ntrt is 3 or 4

  if (!length(S)) {
    S  <- rep(1, length(Y))
    uS <- 1
    nS <- 1
  } else {
    uS <- sort(unique(S))
    nS <- length(uS) 
  }
  Y1 <- Y != 0
 
  # Treatment group 1
  tmp <- getCountVecs(1, trt, Y1, S, uS, nS, op) 
  x1  <- tmp$x
  n1  <- tmp$n

  # Treatment group 2
  tmp <- getCountVecs(2, trt, Y1, S, uS, nS, op) 
  x2  <- tmp$x
  n2  <- tmp$n

  # Control group(s)
  if (ntrt == 3) {
    tmp <- getCountVecs(0, trt, Y1, S, uS, nS, op) 
    x3  <- tmp$x
    n3  <- tmp$n
    x4  <- NA
    n4  <- NA
  } else {
    tmp <- getCountVecs(3, trt, Y1, S, uS, nS, op) 
    x3  <- tmp$x
    n3  <- tmp$n
    tmp <- getCountVecs(4, trt, Y1, S, uS, nS, op) 
    x4  <- tmp$x
    n4  <- tmp$n
  }

  list(x1=x1, x2=2, x3=x3, x4=x4, n1=n1, n2=n2, n3=n3, n4=n4) 
 
} # END: getCounts

checkDataObj <- function(obj, name, vec=1, nrows=0) {

  n <- length(obj)
  if (!n) stop(paste("ERROR: ", name, " has length 0", sep=""))
  if (!is.numeric(obj)) stop(paste("ERROR: ", name, " must be numeric", sep=""))
  vFlag  <- is.vector(obj)
  mFlag  <- is.matrix(obj)
  dfFlag <- is.data.frame(obj)
  if (vFlag + mFlag + dfFlag == 0) stop(paste("ERROR: ", name, " must be a vector, matrix or data frame", sep=""))
  nr <- nrow(obj)
  if (vFlag && !vec) {
    dim(obj) <- c(n, 1)
    nr       <- n
  }
  nc <- ncol(obj)
  if (is.null(nr)) nr <- 0
  if (is.null(nc)) nc <- 0
  if (vec && (nc > 1)) stop(paste("ERROR: ", name, " has an incorrect number of columns", sep=""))
  if (vec) nr <- n
  if (nrows && (nr != nrows)) {
    if (vec) {
      stop(paste("ERROR: ", name, " has the wrong length", sep=""))
    } else {
      stop(paste("ERROR: ", name, " has an incorrect number of rows", sep=""))
    }
  }
  if (dfFlag) {
    tmp <- try(as.matrix(obj), silent=TRUE)
    if ("try-error" %in% class(tmp)) stop(paste("ERROR: ", name, " cannot be coerced to a matrix", sep=""))
  }

  NULL

} # END: checkDataObj

checkVector <- function(obj, name, len=0) {

  n <- length(obj)
  if (!n) stop(paste("ERROR: ", name, " has length 0", sep=""))
  if (!is.numeric(obj)) stop(paste("ERROR: ", name, " must be a numeric vector", sep=""))
  vFlag  <- is.vector(obj)
  if (!vFlag) stop(paste("ERROR: ", name, " must be a numeric vector", sep=""))
  if (len && (len != n)) stop(paste("ERROR: ", name, " has the wrong length", sep=""))
  if (any(!is.finite(obj))) stop(paste("ERROR: ", name, " must have non-negative values", sep=""))
  if (any(obj < 0)) stop(paste("ERROR: ", name, " must have non-negative values", sep=""))

  NULL

} # END: checkVector



checkData <- function(trt, Y, X=NULL, S=NULL) {

  n <- length(Y)
  checkDataObj(Y,   "outcome", vec=1, nrows=0)
  checkDataObj(trt, "trtGroup", vec=1, nrows=n)
  if (length(X))  checkDataObj(X, "covars", vec=0, nrows=n)
  if (length(S))  checkDataObj(S, "strata", vec=1, nrows=n)

  NULL

} # END: checkData


check.options <- function(op, method="counts") {

  valid <- c("alpha", "altParam", "deltaMethod", "LRT", 
             "estRR", "fast", "refPop", "DEBUG", "fast.c")
  op <- default.list(op, valid, 
                     list(0.05, 0, 0, 0,
                          0, 0, NULL, 0, 1))
  refPop <- check.refPop(op[["refPop", exact=TRUE]])
  check.alpha(op$alpha)
  check.binary(op$altParam)
  check.binary(op$deltaMethod)
  check.binary(op$LRT)
  check.binary(op$estRR)
  check.binary(op$fast)
  check.binary(op$fast.c)
  
  nm  <- names(op)
  tmp <- !(nm %in% valid)
  if (any(tmp)) {
    str <- paste(nm[tmp], collapse="," , sep="")
    msg <- paste("ERROR: ", str, " are not valid options", sep="")
    stop(msg)
  } 

  if (method == "adj") {
    if (op$deltaMethod) stop("ERROR: deltaMethod = 1 is not valid")
    if (op$altParam) stop("ERROR: altParam = 1 is not valid")
  } else {
    if (op$estRR) stop("ERROR: estRR = 1 is not valid")
  }

  op$contCorrect <- 0

  op

} # END: check.options

DRRCI_main <- function(x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,
                  tr=NA,Y=NA,Z=NA,alpha=0.05,altParam=0,deltaMethod=0,
                  refPop=NULL,LRT=0,estRR=0,fast=0, op=NULL) {

  DEBUG <- op$DEBUG
  if (!length(refPop) || is.na(refPop)) refPop <- "NULL"
  if (DEBUG) cat(paste("alpha=", alpha, ", altParam=", altParam, ", deltaMethod=", deltaMethod,
             ", refPop=", refPop, ", LRT=", LRT, ", estRR=", estRR, ", fast=", fast, "\n", sep=""))
  if (estRR == 0) {
    out <- DRRCI3(x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                  tr=tr,Y=Y,Z=Z,alpha=alpha,altParam=altParam,
                  deltaMethod=deltaMethod,refPop=refPop,LRT=LRT,
                  estRR=estRR,fast=fast, op=op)
  }
  if (estRR == 1 & is.na(x1[1])) {
    out1 <- DRRCI3(x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,tr=tr,
                   Y=Y,Z=Z,alpha=alpha,altParam=altParam,deltaMethod=deltaMethod,
                   refPop=refPop,LRT=LRT,estRR=estRR, op=op)
    out2 <- DRRCI3(x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                   tr=ifelse(tr==1,2,ifelse(tr==2,1,0)),Y=Y,Z=Z,alpha=alpha,
                   altParam=altParam,deltaMethod=deltaMethod,refPop=refPop,
                   LRT=LRT,estRR=estRR, op=op)
    out <- list(EST=matrix(c(out1$EST,out2$EST),nrow=1),LB=matrix(c(out1$LB,out2$LB),nrow=1),
                UB=matrix(c(out1$UB,out2$UB),nrow=1),beta=out1$beta)
    colnames(out$EST)=colnames(out$LB)=colnames(out$UB)=c("Treatment 2","Treatment 3")
  }
  if (estRR == 1 & !is.na(x1[1])) print("ERROR: DRRCI can only estimate RR when data is entered as Y, Z, and tr")

  out

} # END: DRRCI_main


check.refPop <- function(refPop) {

  m <- length(refPop)
  if (m > 1) stop("ERROR: refPop is not valid")
  if (m) {
    valid  <- c("All", "Treat", "Cont", "Opt", "Ind")
    if (!(refPop %in% valid)) stop("ERROR: refPop is not valid")
  } else {
    refPop <- "NULL"
  }

  refPop

} # END: check.refPop

check.alpha <- function(x) {

  m <- length(x)
  if (!m || (m > 1)) stop("ERROR: alpha is not valid")
  if ((x <= 0) || (x >= 1)) stop("ERROR: alpha must be between 0 and 1")

  NULL

} # END: check.alpha

check.binary <- function(obj, name) {

  n <- length(obj)
  if ((n != 1) || !is.finite(obj) || !(n %in% 0:1)) {
    msg <- paste("ERROR: ", name, " is not valid", sep="")
    stop(msg)
  }
  
  NULL

} # END: check.obj

DRRCI3 <- function(x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,tr=NA,Y=NA,Z=NA,alpha=0.05,
                   altParam=0,deltaMethod=0,refPop=NA,LRT=0,estRR=0,fast=0, op=NULL){

  DEBUG <- op$DEBUG

  if (!is.na(x1[1]) && (LRT==0) && (altParam==0) && (length(x1)>1)  && (refPop == "NULL")){
    x1a=x1
    x2a=x2
    r.simp <- (sum(x2)/sum(n2) - sum(x1)/sum(n1))/(sum(x3)/sum(n3))
    p1.simp <- sum(x1)/sum(n1)
    p2.simp <- sum(x2)/sum(n2)
    if(p1.simp==0) p1.simp= 0.001/sum(n1)
    if(p2.simp==0) p2.simp= 0.001/sum(n2)
    n1 <- ifelse(x1a==0 | x2a==0,n1+0.01/p1.simp,n1)
    n2 <- ifelse(x1a==0 | x2a==0,n2+0.01/p2.simp,n2)
    x1 <- ifelse(x1a==0 | x2a==0,x1+0.01,x1)
    x2 <- ifelse(x1a==0 | x2a==0,x2+0.01,x2)

  }



  ####IDENTIFY METHOD
  unCond=0
  method.rr <- "OP"
  if (!is.na(x4[1]))            method.rr = "TP"
  if (length(Y)>1 & estRR==0)   method.rr = "CO"
  if (length(Y)>1 & estRR==1)   method.rr = "COR"
  if (DEBUG) cat(paste("method.rr=", method.rr, "\n", sep=""))

  ####NECESSARY WARNINGS
  if (!is.na(x1[1])  & length(Y) > 1 ) print("WARNING: You should not provide BOTH a 2x3 table AND disease/treatment/covariates")
  if ((refPop != "NULL") & (length(Y) > 1) ) print("WARNING: You should not provide BOTH a reference population table AND disease/treatment/covariates")
  if (deltaMethod==1 & length(Y) > 1 ) print("WARNING: You should not ask for the delta method AND provide disease/treatment/covariates")
  if (deltaMethod==0){

    init.r <- init.R1 <- NA
    ###Initial values for r
    if (method.rr=="OP") init.r   <-  startR.OP(x1,x2,x3,n1,n2,n3)
    if (method.rr=="TP") init.r   <-  startR.TP(x1,x2,x3,x4,n1,n2,n3,n4)
    if (method.rr=="CO") init.r   <-  startR.CO(tr,Z,Y)
    if (method.rr=="COR")init.R1  <-  startR.COR(tr,Z,Y)


    ####Add intercept to Z if necessary
    #if (length(Y) > 1) {if (length(Z)==length(Y) & sum(Z==1)     != length(Y)) Z <-  cbind(1,Z)
    #if (length(Z)!=length(Y) & sum(Z[,1]==1) != length(Y)) Z <-  cbind(1,Z)
    #if (is.null(colnames(Z))) colnames(Z) <- paste0("Z",c(1:ncol(Z))-1)
    #colnames(Z)[1]="Int"}

    init.pa=init.p3=init.p4=init.t=init.tp=init.k=init.b=init.R2=NA

    if (method.rr != "CO" & method.rr != "COR"){

      ####Create x4 and n4 when needed
      if (sum(!is.na(n4))==0) n4 <- rep(0,length(x1))
      if (sum(!is.na(x4))==0) x4 <- rep(0,length(x1))

      ####Remove bad strata
      ###no individuals in one of the groups
      badStrata1 <- c(1:length(n1))[ifelse(n1==0 | n2==0 | n3==0 | (n4==0 & sum(n4)>0),1,0)==1]
      ###no events in the control group
      badStrata2 <- c(1:length(n1))[x3==0]
      if (altParam==0) badStrata <- sort(unique(c(badStrata1,badStrata2)))
      if (altParam==1) badStrata <- badStrata1
      if (length(badStrata) > 0){
        x1 <- x1[-badStrata]
        x2 <- x2[-badStrata]
        x3 <- x3[-badStrata]
        x4 <- x4[-badStrata]
        n1 <- n1[-badStrata]
        n2 <- n2[-badStrata]
        n3 <- n3[-badStrata]
        n4 <- n4[-badStrata]
      }

      ####Identify reference group if needed
      if (refPop=="All")     nR <- n1+n2+n3+n4
      if (refPop=="Treat")     nR <- n1+n2
      if (refPop=="Cont")     nR <- n3+n4
      if ((refPop=="Opt") & (method.rr=="OP"))   {
        p1.hat  <- sum(x1)/sum(n1)
        p2.hat  <- sum(x2)/sum(n2)
        p3.hat  <- sum(x3)/sum(n3)
        pa.hat  <- (p1.hat+p2.hat)/2
        pie.num <- (pa.hat*(1/n1 + 1/n2) + (0.01+(p1.hat-p2.hat)^2)/(n3*p3.hat))^(-1)
        pie     <- pie.num/sum(pie.num)
        nR      <- pie*sum(n1+n2+n3) 
      }
      if ((refPop=="Opt") && (method.rr=="TP"))   {
        cat("WARNING: The option refPop=Opt does not work for method.rr=TP. Defaults to refPop=All\n")
        nR <- n1+n2+n3+n4
      }


      ###Recalibrates x's and n's to reference group if needed
      if (refPop != "NULL"){
        snR <- sum(nR)
        E1 <- sum((nR/n1)*x1)
        E2 <- sum((nR/n2)*x2)
        E3 <- sum((nR/n3)*x3)
        E4 <- sum((nR/n4)*x4)
        n1 <- (snR^2)/sum(nR^2/n1)
        n2 <- (snR^2)/sum(nR^2/n2)
        n3 <- (snR^2)/sum(nR^2/n3)
        n4 <- (snR^2)/sum(nR^2/n4)
        x1 <- n1*E1/snR
        x2 <- n2*E2/snR
        x3 <- n3*E3/snR
        x4 <- n4*E4/snR
      }

      ###number of strata
      s <- length(x1)
      #goodStrata <- c(1:s)[(x1+x2+x3+x4>0) & (x1+x2 > 0 | altParam==1 | method.rr!="OP")]
      goodStrata <- c(1:s)[(x1+x2+x3+x4>0)]
      if (length(goodStrata) < s & s > 1) warning("WARNING: Dropping strata without any events (or events in treatment)")

      if ((s > 1)  && (length(goodStrata) > 1)) {
        x1 <- x1[goodStrata]
        x2 <- x2[goodStrata]
        x3 <- x3[goodStrata]
        x4 <- x4[goodStrata]
        n1 <- n1[goodStrata]
        n2 <- n2[goodStrata]
        n3 <- n3[goodStrata]
        n4 <- n4[goodStrata]
      }
      s <- length(x1)


      N       <- n1+n2+n3+n4
      init.p3 <-  x3/n3
      init.p3 <- ifelse(init.p3 < 0.000001,0.000001,ifelse(init.p3 > 0.9999,0.9999,init.p3))
      init.p4 <-  x4/n4
      init.tp <-  (x1/n1)/(x3/n3)+(x2/n2)/(x4/n4)
      init.pa <- (x1/n1+x2/n2)/2
      init.pa <- ifelse(init.pa < 0.000001,0.000001,ifelse(init.pa > 0.9999,0.9999,init.pa))
      init.pa <- ifelse(init.pa < init.r*init.p3/2+0.000001,init.r*init.p3/2+0.000001,init.pa)
      init.pa <- ifelse(init.pa < -init.r*init.p3/2+0.000001,-init.r*init.p3/2+0.000001,init.pa)

      init.t  <- sum(N*init.pa)/sum(N*init.p3)
      tr <- Y <- Z <- NA
    }

    if (substr(method.rr,1,2) == "CO"){
      ncz <- ncol(Z)
      if (ncz > 1) {
        badCol <- rep(0,ncz)
        for (i in 2:ncz) if (sum(Y[abs(Z[,i]) > 0])==0) badCol[i] <- 1
        Z <- Z[,badCol==0, drop=FALSE]
      }
      init.b2 <- log(mean(Y[tr==0],na.rm=T))
      init.b  <- c(ifelse(is.na(init.b2) | init.b2 < -15, -15,init.b2),rep(0,ncol(Z)-1))
      init.k  <- mean(Y[tr!=0],na.rm=T)/exp(init.b[1])
      if (method.rr == "CO") init.k  <- ifelse(init.k-init.r/2 < 0.000001, init.r/2 + 0.000001, init.k )
      init.R2 <- mean(Y[tr==2],na.rm=T)/mean(Y[tr==0],na.rm=T)
      x1 <- x2 <- x3 <- x4 <- n1 <- n2 <- n3 <- n4 <- NA
    }

    EST <- LB <- UB <- NA
    if ( ((sum(x3, na.rm=TRUE)>=0) && ((sum(x4, na.rm=TRUE)>=0) || (sum(n4, na.rm=TRUE)==0))) || (substr(method.rr,1,2) == "CO") ){
      if (method.rr=="OP") {
        if (DEBUG) cat("calling getEst.OP\n")
        EST <- getEst.OP(init.pa,init.t,init.p3,init.r,x1,x2,x3,n1,n2,n3,altParam,unCond, op=op)
      } else if (method.rr=="TP") {
        if (DEBUG) cat("calling getEst.TP\n")
        EST <- getEst.TP(init.p4,init.tp,init.p3,init.r,x1,x2,x3,x4,n1,n2,n3,n4,op=op)
      } else if (method.rr=="CO") {
        if (DEBUG) cat("calling getEst.CO\n")
        EST <- getEst.CO(init.k,init.b,init.r,Y,Z,tr, op=op)
      } else if (method.rr=="COR") {
        if (DEBUG) cat("calling getEst.COR\n")
        EST <- getEst.COR(init.R2,init.b,init.R1,Y,Z,tr, op=op)
      }
    }

    if ( ((sum(x3, na.rm=TRUE)==0) || ((sum(x4, na.rm=TRUE)==0) && (sum(n4, na.rm=TRUE)>0))) && (substr(method.rr,1,2) != "CO"))  EST     <- 0
    if (!is.na(EST) && (LRT==0)) {
      if (DEBUG) cat("Calling boundR\n")
#print("LOWER BOUND")

      LB <- boundR(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,alpha2=alpha/2,type="LB",altParam=altParam,
                   method.rr=method.rr,unCond=unCond,op=op)
      if (DEBUG) cat("Calling boundR\n")
#print("UPPER BOUND")
      UB <- boundR(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,alpha2=alpha/2,type="UB",altParam=altParam,
                   method.rr=method.rr,unCond=unCond,op=op)
    }
    if (!is.na(EST) && (LRT==1)) {
      if (DEBUG) cat("Calling boundR.LRT\n")
#print("LOWER BOUND")
      try(LB <- boundR.LRT(EST,x1,x2,x3,x4,n1,n2,n3,n4,Y,Z,tr,init.pa,init.p3,init.p4,init.t,init.tp,init.k,
                           init.b,init.R1,init.R2,method.rr=method.rr,altParam=altParam,alpha2=alpha/2,
                           type="LB", op=op), silent=TRUE)
#print("UPPER BOUND")

      if (DEBUG) cat("Calling boundR.LRT\n")
      try(UB <- boundR.LRT(EST,x1,x2,x3,x4,n1,n2,n3,n4,Y,Z,tr,init.pa,init.p3,init.p4,init.t,init.tp,init.k,
                           init.b,init.R1,init.R2,method.rr=method.rr,altParam=altParam,alpha2=alpha/2,
                           type="UB", op=op), silent=TRUE)
    }

    beta               <- NULL
    if (substr(method.rr,1,2) == "CO"){
      cname              <- colnames(Z)
      beta               <- matrix(ncol=3, nrow=length(cname))
      if (length(cname)==0) cname <- rep("--",ncol(Z))
      bestEst            <- NA
      if (method.rr == "CO")  try(bestEst         <- getEst.CO( init.k,init.b,init.r,Y,Z,tr,allRes=1, op=op), silent=TRUE)
      if (method.rr == "COR") try(bestEst         <- getEst.COR(init.R2,init.b,init.R1,Y,Z,tr,allRes=1, op=op), silent=TRUE)
      if (length(bestEst) > 1){
        lbe               <- length(bestEst)
        if (method.rr == "CO")  III               <- fI.CO(bestEst[lbe], bestEst[1],bestEst[-c(1,lbe)],tr,length(Y),Z)
        if (method.rr == "COR") III               <- fI.COR(bestEst[lbe],bestEst[1],bestEst[-c(1,lbe)],tr,length(Y),Z)
        SE                <- rep(NA,lbe-2)
        try(SE            <- sqrt(diag(solve(III)))[-c(1,lbe)], silent=TRUE)
        betas             <- bestEst[-c(1,lbe)]
        beta              <- data.frame("name"=cname,"betas"=betas,"SE"=SE, stringsAsFactors=FALSE)
        colnames(beta)    <- c("","Beta","SE")
      }
    }

    out <- list(EST=EST,LB=LB,UB=UB)
    if (length(beta)) out$beta <- beta

    ####end if deltaMethod==0
  }


  if (deltaMethod==1) out <- deltaMethod(x1=x1,x2=x2,x3=x3,n1=n1,n2=n2,n3=n3,alpha=alpha,refPop=refPop)

  out

}

optim.force.C <-  function(pa=NA,p3=NA,t=NA,r=NA,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="OP",
                         x1,x2,x3,x4,n1,n2,n3,n4,n,Z,Y,tr,fix.r=0,allRes=0, DEBUG=0, 
                         initStepSize=1.0, stopTol=1e-3){

  if (DEBUG) cat("Begin: optim.force.C\n")
  debug  <- DEBUG
  s      <- length(x1)
  OP0  <- (method.rr == "OP") && !altParam
  OP1  <- (method.rr == "OP") && altParam
  OP   <- OP0 || OP1
  CO   <- method.rr == "CO"   
  COR  <- method.rr == "COR"
  TP   <- method.rr == "TP"

  # Check objects pa, p3 etc  Sometime they do not exist
  f0 <- -9999999
  if (OP0) {
    allParam  <- c(pa,p3,r)
    paramType <- c(rep("pa",length(pa)),rep("p3",length(p3)),rep("r",length(r)))
    method    <- 0
  } else if (OP1) {
    allParam  <- c(t,p3,r)
    paramType <- c(rep("t",length(t)),  rep("p3",length(p3)),rep("r",length(r)))
    method    <- 1
  } else if (CO) {
    allParam  <- c(k,b,r)
    paramType <- c(rep("k",length(k)),  rep("b",length(b)),  rep("r",length(r)))
    method    <- 2
  } else if (COR) {
    allParam  <- c(R2,b,R1)
    paramType <- c(rep("R2",length(R2)),rep("b",length(b)),  rep("R1",length(R1)))
    method    <- 3
  } else if (TP) {
    allParam  <- c(tp,p3,p4,r)
    paramType <- c(rep("tp",length(tp)),rep("p3",length(p3)),rep("p4",length(p4)),rep("r",length(r)))
    method    <- 4
  }  

  nallParam    <- length(allParam)
  out          <- NA
  vec1         <- c(1, nallParam)
  maxiter      <- 10000
  probEps      <- 1.0e-6
  npossVals    <- 10
  retParms     <- rep(-9999, nallParam)
  retOK        <- -1
  retMin       <- 1e100
  parmOrder    <- 0:(nallParam - 1) # adjusted for C
  ptype        <- gsub("p", "", paramType, fixed=TRUE)
  ptype        <- gsub("R", "", ptype, fixed=TRUE)
  ptype        <- paste(ptype, collapse="", sep="")
  minParmValue <- rep(-9999, nallParam)
  maxParmValue <- rep(-9999, nallParam)
  
  #tmp          <- getMinMaxValues(paramType, probEps)
  #minParmValue <- tmp$minv
  #maxParmValue <- tmp$maxv

  if (OP0 || OP1 || TP) {
    n1x1  <- n1 - x1
    n2x2  <- n2 - x2
    n3x3  <- n3 - x3
  } else {
    x1 <- x2 <- x3 <- n1x1 <- n2x2 <- n3x3 <- 0
  }
  if (TP) {
    n4x4 <- n4 - x4
  } else {
    x4 <- n4x4 <- 0
  }
  if (CO || COR) {
    nsubs   <- length(Y)
    ncovars <- ncol(Z) - 1
  } else {
    Y <- tr <- Z <- nsubs <- ncovars <- 0
  }
  if (!ncovars) Z <- 0

  if (debug > 1) {
    print(paramType)
    print(cbind(allParam, minParmValue, maxParmValue))
    print(paste("f0 = ", f0, ", fix.r = ", fix.r, sep=""))
  }

  # Integer valued scalars
  iargs <- c(ncovars, nsubs, maxiter, method, fix.r, s, nallParam, npossVals, debug)
  dargs <- c(probEps, initStepSize, stopTol, f0)
  tmp <- try(.C("C_optim_force", as.integer(iargs), as.numeric(dargs), 
            as.numeric(x1), as.numeric(x2), as.numeric(x3), as.numeric(x4), 
            as.numeric(n1x1), as.numeric(n2x2), as.numeric(n3x3), as.numeric(n4x4),
            as.integer(Y), as.integer(tr), as.numeric(t(Z)), 
            as.numeric(allParam), as.numeric(minParmValue), as.numeric(maxParmValue), 
            as.integer(parmOrder), as.character(ptype),
            retParms=as.numeric(retParms), retOK=as.integer(retOK), 
            #retMin=as.numeric(retMin)), silent=TRUE)
            retMin=as.numeric(retMin), PACKAGE="DiffRelRisk"), silent=TRUE)
  if ("try-error" %in% class(tmp)) {
    if (debug) print(tmp)
    stop("ERROR in C_optim_force")
  }

  retOK    <- tmp$retOK
  retParms <- tmp$retParms
  retMin   <- tmp$retMin
  allParam <- retParms
 
  if (retOK){

    if (OP || TP) {
      oneTos     <- 1:s
      sp1To2s    <- (s+1):(2*s)
      twosp1     <- 2*s+1
      sp2        <- s + 2
      twoTosp1   <- 2:(s+1)
      twosp1To3s <- (2*s+1):(3*s)
    }

    if (OP0) {
      if (!allRes) {
        if (!fix.r) {    
          out <- allParam[twosp1]
        } else {
          out <- list("pa"=allParam[oneTos],"p3"=allParam[sp1To2s],"r"=allParam[twosp1],"t"=NA)
        }
      } else {
        if (fix.r) {
          out <- c(allParam[oneTos],allParam[sp1To2s])
        } else {
          out <- c(allParam[oneTos],allParam[sp1To2s],allParam[twosp1])
        }
      }
    } else if (OP1) {
      if (!allRes) {
        if (!fix.r) {    
          out <- allParam[sp2]
        } else {
          out <- list("pa"=NA, "p3"=allParam[twoTosp1], "r"=allParam[sp2], "t"=allParam[1])
        }
      } else {
        if (fix.r) {
          out <- c(allParam[1],allParam[twoTosp1])
        } else {
          out <- c(allParam[1],allParam[twoTosp1],allParam[sp2])
        }
      }
    } else if (CO) {
      if (!allRes) {
        if (!fix.r) {    
          out <- allParam[nallParam]
        } else {
          out <- list("r"=allParam[nallParam],"k"=allParam[1],"b"=allParam[-vec1])
        }
      } else {
        if (fix.r) {
          out <- allParam[-nallParam]
        } else {
          out <- allParam
        }
      }
    } else if (COR) {
      if (!allRes) {
        if (!fix.r) {    
          out <- allParam[nallParam]
        } else {
          out <- list("R1"=allParam[nallParam],"R2"=allParam[1],"b"=allParam[-vec1])
        }
      } else {
        if (fix.r) {
          out <- allParam[-nallParam]
        } else {
          out <- allParam
        }
      }
    } else if (TP) {
      if (!allRes) {
        if (!fix.r) {    
          out <- allParam[nallParam]
        } else {                        
          out <- list("t"=allParam[oneTos],"r"=allParam[nallParam],"p01"=allParam[sp1To2s],"p02"=allParam[twosp1To3s])
        }
      } else {
        if (fix.r) {
          out <- allParam[-nallParam]
        } else {
          out <- allParam
        }
      }
    }
  } # END: if (mult < 0.1)

  ret <- list("out"=out,"ok"=retOK, retMin=retMin)

  if (DEBUG > 1) {
    cat("optim.force.C return object:\n")
    print(ret) 
  } 
  if (DEBUG) cat("End: optim.force.C\n")

  ret

} # END: optim.force.C



####################################################################################################################
#####
#####Delta Method
#####
####################################################################################################################

deltaMethodCheckCnt <- function(v, newValue=0.5) {

  flag <- FALSE
  tmp  <- v == 0
  if (any(tmp)) {
    v[tmp] <- newValue
    flag   <- 1
  } 

  list(flag=flag, counts=v)

} # END: deltaMethodCheckCnt

deltaMethodCheckCounts <- function(x1, x2, x3, n1, n2, n3) {

  changed <- FALSE
  tmp     <- deltaMethodCheckCnt(x1)
  x1      <- tmp$counts
  changed <- changed || tmp$flag
  tmp     <- deltaMethodCheckCnt(x2)
  x2      <- tmp$counts
  changed <- changed || tmp$flag
  tmp     <- deltaMethodCheckCnt(x3)
  x3      <- tmp$counts
  changed <- changed || tmp$flag
  tmp     <- deltaMethodCheckCnt(n1)
  n1      <- tmp$counts
  changed <- changed || tmp$flag
  tmp     <- deltaMethodCheckCnt(n2)
  n2      <- tmp$counts
  changed <- changed || tmp$flag
  tmp     <- deltaMethodCheckCnt(n3)
  n3      <- tmp$counts
  changed <- changed || tmp$flag

  list(changed=changed, x1=x1, x2=x2, x3=x3, n1=n1, n2=n2, n3=n3)

} # END: deltaMethodCheckCounts

deltaMethod <- function(x1=NA,x2=NA,x3=NA,n1=NA,n2=NA,n3=NA,alpha=0.05,refPop=NA){

  # For the variances in the delta-method, you can add 0.5 to each group w/ 0.

  ret <- deltaMethodF2(x1=x1,x2=x2,x3=x3,n1=n1,n2=n2,n3=n3,alpha=alpha,refPop=refPop)
  
  # Check the counts
  tmp <- deltaMethodCheckCounts(x1, x2, x3, n1, n2, n3)
  if (tmp$changed) {
    ret2   <- deltaMethodF2(x1=tmp$x1,x2=tmp$x2,x3=tmp$x3,n1=tmp$n1,n2=tmp$n2,n3=tmp$n3,
                           alpha=alpha,refPop=refPop)
    ret$LB <- ret2$LB
    ret$UB <- ret2$UB
  }
 
  ret

} # END: deltaMethod

deltaMethodF2 <- function(x1=NA,x2=NA,x3=NA,n1=NA,n2=NA,n3=NA,alpha=0.05,refPop=NA){


  beta <- NULL
  za   <- qnorm(1-alpha/2)

  if (length(x1)==1){
    p1.hat <- x1/n1
    p2.hat <- x2/n2
    p3.hat <- x3/n3

    s1.hat <- p1.hat*(1-p1.hat)/n1
    s2.hat <- p2.hat*(1-p2.hat)/n2
    s3.hat <- p3.hat*(1-p3.hat)/n3


    EST  <- c((p2.hat-p1.hat)/p3.hat)
    SE   <- sqrt(s1.hat/p3.hat^2 + s2.hat/p3.hat^2 + (p2.hat-p1.hat)^2*s3.hat/p3.hat^4)
    LB   <- c(EST-za*SE)
    UB   <- c(EST+za*SE)
  }

  if ( (length(x1)>1) && (refPop != "Ind") ){

    if (refPop=="NULL")   refPop <- "All"
    if (refPop=="All")   nR <- n1+n2+n3
    if (refPop=="Treat") nR <- n1+n2
    if (refPop=="Cont")  nR <- n3
    if (refPop=="Opt")   {
      p1.hat  <- sum(x1)/sum(n1)
      p2.hat  <- sum(x2)/sum(n2)
      p3.hat  <- sum(x3)/sum(n3)
      pa.hat  <- (p1.hat+p2.hat)/2
      pie.num <- (pa.hat*(1/n1 + 1/n2) + (0.01+(p1.hat-p2.hat)^2)/(n3*p3.hat))^(-1)
      pie     <- pie.num/sum(pie.num)
      nR      <- pie*sum(n1+n2+n3) }

    p1.hat <- x1/n1
    p2.hat <- x2/n2
    p3.hat <- x3/n3

    pie.m <- nR/sum(nR)
    p.star.1 <- sum(pie.m*p1.hat)
    p.star.2 <- sum(pie.m*p2.hat)
    p.star.3 <- sum(pie.m*p3.hat)


    s1.hat <- sum(pie.m^2*p1.hat*(1-p1.hat)/n1)
    s2.hat <- sum(pie.m^2*p2.hat*(1-p2.hat)/n2)
    s3.hat <- sum(pie.m^2*p3.hat*(1-p3.hat)/n3)

    EST <- (p.star.2-p.star.1)/p.star.3
    SE   <- sqrt(s1.hat/p.star.3^2 + s2.hat/p.star.3^2 + s3.hat*(p.star.2-p.star.1)^2/p.star.3^4)
    LB   <- c(EST-za*SE)
    UB   <- c(EST+za*SE)
  }


  if ((length(x1)>1) && (refPop == "Ind")){


    p1.hat <- x1/n1
    p2.hat <- x2/n2
    p3.hat <- x3/n3
    x01.t <- sum(n1*p3.hat)
    x02.t <- sum(n2*p3.hat)

    s1.hat <- sum(n1*p1.hat*(1-p1.hat))
    s2.hat <- sum(n2*p2.hat*(1-p2.hat))

    s01.hat  <- sum(n1^2*p3.hat*(1-p3.hat)/n3)
    s02.hat  <- sum(n2^2*p3.hat*(1-p3.hat)/n3)
    s012.hat <- sum(n2*n1*p3.hat*(1-p3.hat)/n3)

    EST <- sum(x2)/sum(x02.t) - sum(x1)/sum(x01.t)
    SE   <- sqrt(s1.hat/x01.t^2 + s2.hat/x02.t^2 + sum(x1)^2*s01.hat/x01.t^4 + sum(x2)^2*s02.hat/x02.t^4 -2*sum(x2)*sum(x1)*s012.hat/(x01.t^2*x02.t^2))
    LB   <- c(EST-za*SE)
    UB   <- c(EST+za*SE)
  }

  list(EST=EST,LB=LB,UB=UB)

} # END: deltaMethodF2

call_findR.int <- function(r=NA,alpha2=0.025,x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,
                                tr=NA,Z=NA,Y=NA,altParam=0,method.rr="OP",unCond=0,
                                 sv=NA,use.force=1, op=NULL) {
    BADVALUE <- 99999

    out2 <- try(findR.int(r=r,alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                                tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,
                                 sv=sv,use.force=use.force, op=op),silent=TRUE)
    flag <- ("try-error" %in% class(out2)) || !is.finite(out2$out) || (abs(out2$out) > BADVALUE)

    if (flag && !all(is.na(sv))) {
      out2 <- try(findR.int(r=r,alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                                tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,
                                 sv=NA,use.force=use.force, op=op),silent=TRUE)
#print("call_findR.int")
#print(out2)
      flag <- ("try-error" %in% class(out2)) || !is.finite(out2$out) || (abs(out2$out) > BADVALUE)
    }

    if (flag) {
      return("findR.int failed")
    }

    out2

} # END: call_findR.int

getBounds.check <- function(list1, list2, type, mult, alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                                tr=tr,Z=Z,Y=Y,altParam=0,method.rr="OP",unCond=0,sv=NA, op=NULL) {

  # Function to check the new return object after initial call that bounds the point we want
  # Sometimes when subdividing the interval, the derivative will slighly decrease then increase
  #  in a region that does not contain the point we want. When this happens, we will go back to
  #  the midpoint of the original interval and start from there.

  der1 <- list1$POSS.D
  if (!is.list(list2)) list2 <- list(POSS.D=rep(1e100, 3))
  der2 <- list2[["POSS.D", exact=TRUE]]
  if (!length(der2)) der2 <- rep(1e100, 3)
  der1.2 <- der1[2] 
  if (all(is.finite(der2)) && (min(der2) < der1.2)) return(list2)

  sv       <- (list1$SVLIST)[[1]]
  pt1      <- list1$POSS.R
  mid      <- pt1[2]
  N        <- 5
  h        <- abs(pt1[3] - mid)/N

#print(paste("111  h=", h, sep=""))
#print(cbind(list1$POSS.R, list1$POSS.D, list2$POSS.R, list2$POSS.D))



  for (i in 1:(N-1)) {
    mult  <- i*h
    min.r <- mid - mult - mult  # Because getBounds.force will add mult to min.r for the first point
    max.r <- mid + mult 
#print(paste("min.r=", min.r, ", max.r=", max.r, sep=""))
    tmp   <- getBounds.force(min.r, max.r, "UB", mult, alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                             tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,sv=sv, op=op)
    if (!is.list(tmp)) {     
      next
    }
    POSS.D <- tmp$POSS.D
    if (any(POSS.D <= der1.2)) return(tmp)
  }

#print("222")

  return("getBounds.check failed")

} # END: getBounds.check

getBounds.force <- function(min.r, max.r, type, mult, alpha2=0.025,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                                tr=tr,Z=Z,Y=Y,altParam=0,method.rr="OP",unCond=0,sv=NA, op=NULL) {
  POSS.D   <- rep(NA, 3)
  POSS.R   <- rep(NA, 3)
  SVLIST   <- list(NA, NA, NA)
  UB       <- type == "UB"
  BADVALUE <- 99999
  aDec     <- 0
  aInc     <- 0
  okcnt    <- 0

  # Initialize the first point and change step size if needed
  # We start near the MLE and go to more extreme values
  if (UB) {
    poss.r <- min.r
    #if (poss.r + mult > max.r) mult <- (max.r - poss.r)/10
  } else {
    poss.r <- max.r
    #if (poss.r - mult < min.r) mult <- (poss.r - min.r)/10
  }

  while(1) {
    if (UB) {
      poss.r <- poss.r + mult
      if (poss.r > max.r) return("poss.r > max.r")  
    } else {
      poss.r <- poss.r - mult
      if (poss.r < min.r) return("poss.r < min.r")
    }

    out2 <- call_findR.int(r=poss.r,alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                                tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,
                                 sv=sv,use.force=2, op=op)
    if (!is.list(out2)) {
      return("findR.int failed")
    }


    value <- out2$out

    if (!is.finite(value) || (abs(value) > BADVALUE)) {
      return("Bad return value from findR.int")
    }

    sv          <- out2$sv2
    POSS.D      <- c(POSS.D[-1], value)
    POSS.R      <- c(POSS.R[-1], poss.r)
    SVLIST      <- SVLIST[-1]
    SVLIST[[3]] <- sv
    okcnt       <- okcnt + 1
    possd2      <- POSS.D[2]
   
    if (okcnt > 1) {
#print(POSS.R)
#print(POSS.D)
      if (value < possd2) {
        aDec <- 1
      } else if (value > possd2) {
        if (aDec && (okcnt > 2)) {
          sv <- SVLIST[[1]]
          break 
        } else if (okcnt > 100) {
          return("not decreasing")
        }
      } else {
        return("new value == prev value")
      }
    }
  }

  list(POSS.R=POSS.R, POSS.D=POSS.D, sv=sv, SVLIST=SVLIST)

} # END: getBounds.force

boundR2.force <- function(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,R2=NA,alpha2=0.025,type="UB",
                       altParam=0,method.rr="OP",unCond=0, op=NULL) {

  DEBUG     <- op$DEBUG
  out       <- NA
  if (!op$fast) {
    multVec <- c(0.1, 0.05, 0.01, 0.005)
    if (abs(EST) <= 0.05) multVec <- c(multVec, 0.001)
  } else {
    multVec <- c(0.05, 0.01)
  }
  NTRY      <- length(multVec)
  BADVALUE  <- 99999
  CONST     <- 3
  ntry      <- 0
  eps       <- 1e-3
  
  if (type=="UB") {
    UB      <- 1
    min.r   <- EST + eps
    max.r   <- EST + CONST
  } else {
    UB      <- 0
    min.r   <- EST - CONST
    max.r   <- EST - eps
    if (method.rr=="COR") {
      min.r <- max(min.r, eps)
      max.r <- max(max.r, eps)
    }
  }
  minReq2 <- ifelse(!is.na(x1[1]) & (sum(x1==0)>0 | sum(x2==0)>0) ,3, ifelse(is.na(x1[1])  & (sum(Y[tr==1])==0 | sum(Y[tr==2])==0),3,1))
  minReq  <- ifelse(!is.na(x1[1]) & (sum(x1==0)>0 | sum(x2==0)>0) ,1, ifelse(is.na(x1[1])  & (sum(Y[tr==1])==0 | sum(Y[tr==2])==0),3,0.25))
  sv      <- NA
  if (DEBUG) print(paste("min.r=", min.r, ", max.r=", max.r, ", minReq2=", minReq2, ", minReq=", minReq, sep=""))
  while ( (ntry <= NTRY-1) && is.na(out)) {
    ntry <- ntry+1
    sv   <- NA 
    mult <- multVec[ntry]
    tmp  <-getBounds.force(min.r, max.r, type, mult, alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                            tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,sv=sv, op=op)
    if (DEBUG) {
      print(paste("1: try=", ntry, ", mult=", mult, ", getBounds.force return object:", sep=""))
      print(tmp)
    }

    if (!is.list(tmp)) next
    POSS.R <- tmp$POSS.R
    SVLIST <- tmp$SVLIST
#print("1111")
#print(POSS.R)
#print(tmp$POSS.D)
    if (UB) {
      sv <- SVLIST[[1]]
    } else {
      sv <- SVLIST[[1]]
    }
    mult <- mult*0.1
    tmp0 <- tmp
    tmp  <- getBounds.force(min(POSS.R), max(POSS.R), type, mult, alpha2=alpha2,
                            x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                            tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,
                            sv=sv, op=op)
    if (DEBUG) {
      print(paste("2: try=", ntry, ", mult=", mult, ", getBounds.force return object:", sep=""))
      print(tmp)
    }

    tmp  <- getBounds.check(tmp0, tmp, type, mult, alpha2=alpha2,
                            x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                            tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,
                            sv=sv, op=op)

    if (!is.list(tmp)) {
#print(tmp)
      next
    }
    POSS.R <- tmp$POSS.R
    SVLIST <- tmp$SVLIST
#print("2222")
#print(POSS.R)
#print(tmp$POSS.D)

    if (UB) {
      sv <- SVLIST[[1]]
    } else {
      sv <- SVLIST[[1]]
    }

    mult <- mult*0.01
    tmp  <- getBounds.force(min(POSS.R), max(POSS.R), type, mult, alpha2=alpha2,
                            x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,
                            tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,
                            sv=sv, op=op)
    if (DEBUG) {
      print(paste("3: try=", ntry, ", mult=", mult, ", getBounds.force return object:", sep=""))
      print(tmp)
    }

    if (!is.list(tmp)) {
#print(tmp)
      next
    }
    POSS.R <- tmp$POSS.R
    POSS.D <- tmp$POSS.D

#print("3333")
#print(POSS.R)
#print(POSS.D)


    out    <- ifelse(POSS.D[2] > minReq, NA, POSS.R[2])
  }

  out

} # END: boundR2.force





#######################################################################################################################################################################
###
#######################################################################################################################################################################
# Function to create a vector from a matrix
makeVector <- function(x) {

  d <- dim(x)
  if (is.null(d)) return(x)

  nn <- NULL
  if (d[1] == 1) {
    nn <- colnames(x)
  } else if (d[2] == 1) {
    nn <- rownames(x)
  }
  dim(x) <- NULL
  if (!is.null(nn)) names(x) <- nn
  if ((!is.vector(x)) && (is.list(x))) x <- unlist(x)
  x
 
} # END: makeVector

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list



