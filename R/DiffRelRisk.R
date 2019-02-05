#' Difference in Relative Risks
#'
#' This function estimates and provides a confidence interval for the difference between
#' two relative risks. The first option is there are three populations and the goal to estimate
#' the difference: p2/p3 - p1/p3, where p1, p2, and p3 are the risks in each population.
#' The second option is there are four populations and the goal is to estimate the difference:
#' p2/p4 - p1/p3.
#'
#' @import nleqslv
#' @param x1 The number of events in the first group  (either a single number if there is one strata or a vector if there are multiple strata)
#' @param x2 The number of events in the second group (either a single number if there is one strata or a vector if there are multiple strata)
#' @param x3 The number of events in the third group (either a single number if there is one strata or a vector if there are multiple strata)
#' @param x4 The number of events in the fourth group if the study included four groups (either missing if a three-group study, a single number if there is one strata or a vector if there are multiple strata; note: not all options are currently available when there are four groups)
#' @param n1 The total number of subjects in the first group  (either a single number if there is one strata or a vector if there are multiple strata)
#' @param n2 The total number of subjects in the second group (either a single number if there is one strata or a vector if there are multiple strata)
#' @param n3 The total number of subjects in the third group (either a single number if there is one strata or a vector if there are multiple strata)
#' @param n4 The total number of subjects in the fourth group if the study included four groups (either missing if a three-group study, a single number if there is one strata or a vector if there are multiple strata; note: not all options are currently available when there are four groups)
#' @param alpha The two-sided error rate for the confidence interval (i.e. 0.05 for a 95\% confidence interval)
#' @param altParam A binary variable indicating whether to use the full likelihood (altParam=0) or assume the log-linear model holds true (altParam==1) (Note: altParam=1 is equivalent to entering subject Y, tr, Z w/ Z set to binary indicators for strata)
#' @param deltaMethod A binary variable indicating whether the delta method should be used (note: if deltaMethod==1, then default refPop="All")
#' @param refPop The target population when standardization is used (options are "All" -- all subjects, "Treat" -- treated population, "Cont" -- control population; "Opt" -- an optimally chosen population; "Ind" -- indirect standardization for the delta method)
#' @param estRR  A binary variable indicating whether to calculate the confidence intervals for each relative risk, as opposed to the difference (set to 0 for calculating the CI for difference, 1 for calculating the CIs for RR; note only works when data is entered as Z, Y, tr)
#' @param LRT A binary variable indicating whether a likelihood-based confidence interval should be used (set to 0 for score-based confidence interval; set to 1 for likelihood-based confidence interval)
#' @param tr The treatment group for each subject (a vector with one value per subject; values are either 0, 1, or 2; only required if want to adjust for covariates)
#' @param Y The outcome for each subject (a vector with one value per subject; values are either 0 or 1; only required if want to adjust for covariates)
#' @param Z The covariates for each subject (a matrix with one value per subject and one column per covariate; only required if want to adjust for covariates)
#' @param fast A binary variable to indicate that the optimizations should be performed over a coarser grid (set to 1 to speed up analysis; use with caution)
#' @return EST The estimated difference in relative risks
#' @return LB The lower bound for the difference
#' @return UB The upper bound for the difference
#' @return beta The MLE (and asymptotic SE) for the coefficients for individual level covariates (when adjusting for covariates)
#' @examples DRRCI(x1=20,x2=50,x3=100,n1=1000,n2=1000,n3=1000)
#' @examples #The simplest case where I have three groups and
#' @examples #want to estimate (50/1000 - 20/1000)/(100/1000)
#' @examples ##
#' @examples DRRCI(x1=c(15,5),x2=c(30,20),x3=c(50,50),n1=c(600,400),n2=c(400,600),n3=c(500,500))
#' @examples #I now have observed individuals in two strata
#' @examples #here I could use altParam = 1 if I wanted to assume that log-linear model held
#' @examples #here I could use refPop   = "All" if I wanted to use direct standardization
#' @examples ##
#' @examples DRRCI(x1=20,x2=50,x3=100,x4=90,n1=1000,n2=1000,n3=1000,n4=90)
#' @examples #I now have a separate "control" group for each treatment group
#' @examples ##
#' @examples tr = rep(c(0,1,2),each=100)
#' @examples Y  = c(rep(1,30),rep(0,70),rep(1,10),rep(0,90),rep(1,20),rep(0,80))
#' @examples Z  = cbind(1,rnorm(300),rnorm(300))
#' @examples DRRCI(tr=tr,Y=Y,Z=Z)
#' @examples #I now adjust for covariates
#' @export
#'


DRRCI <- function(x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,tr=NA,Y=NA,Z=NA,alpha=0.05,altParam=0,deltaMethod=0,refPop=NA,LRT=0,estRR=0,fast=0){
  if (estRR == 0) out <- DRRCI2(x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,tr=tr,Y=Y,Z=Z,alpha=alpha,altParam=altParam,deltaMethod=deltaMethod,refPop=refPop,LRT=LRT,estRR=estRR,fast=fast)
  if (estRR == 1 & is.na(x1[1])) {out1 <- DRRCI2(x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,tr=tr,Y=Y,Z=Z,alpha=alpha,altParam=altParam,deltaMethod=deltaMethod,refPop=refPop,LRT=LRT,estRR=estRR)
  out2 <- DRRCI2(x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,tr=ifelse(tr==1,2,ifelse(tr==2,1,0)),Y=Y,Z=Z,alpha=alpha,altParam=altParam,deltaMethod=deltaMethod,refPop=refPop,LRT=LRT,estRR=estRR)
  out <- list("EST"=matrix(c(out1$EST,out2$EST),nrow=1),"LB"=matrix(c(out1$LB,out2$LB),nrow=1),"UB"=matrix(c(out1$UB,out2$UB),nrow=1),"beta"=out1$beta)
  colnames(out$EST)=colnames(out$LB)=colnames(out$UB)=c("Treatment 1","Treatment 2")}
  if (estRR == 1 & !is.na(x1[1])) print("ERROR: DRRCI can only estimate RR when data is entered as Y, Z, and tr")
  out
}




DRRCI2 <- function(x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,tr=NA,Y=NA,Z=NA,alpha=0.05,altParam=0,deltaMethod=0,refPop=NA,LRT=0,estRR=0,fast=0){


  if (!is.na(x1[1]) & is.na(refPop) & LRT==0 & altParam==0 & length(x1)>1 & refPop != "Ind"){
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

  ####NECESSARY WARNINGS
  if (!is.na(x1[1])  & length(Y) > 1 ) print("WARNING: You should not provide BOTH a 2x3 table AND disease/treatment/covariates")
  if (!is.na(refPop) & length(Y) > 1 ) print("WARNING: You should not provide BOTH a reference population table AND disease/treatment/covariates")
  if (deltaMethod==1 & length(Y) > 1 ) print("WARNING: You should not ask for the delta method AND provide disease/treatment/covariates")
  if (deltaMethod==0){

    init.r <- init.R1 <- NA
    ###Initial values for r
    if (method.rr=="OP") init.r   <-  startR.OP(x1,x2,x3,n1,n2,n3)
    if (method.rr=="TP") init.r   <-  startR.TP(x1,x2,x3,x4,n1,n2,n3,n4)
    if (method.rr=="CO") init.r   <-  startR.CO(tr,Z,Y)
    if (method.rr=="COR")init.R1  <-  startR.COR(tr,Z,Y)


    ####Add intercept to Z if necessary
    if (length(Y) > 1) {if (length(Z)==length(Y) & sum(Z==1)     != length(Y)) Z <-  cbind(1,Z)
    if (length(Z)!=length(Y) & sum(Z[,1]==1) != length(Y)) Z <-  cbind(1,Z)
    if (is.null(colnames(Z))) colnames(Z) <- paste0("Z",c(1:ncol(Z))-1)
    colnames(Z)[1]="Int"}

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
      if (refPop=="All"      & !is.na(refPop))     nR <- n1+n2+n3+n4
      if (refPop=="Treat"    & !is.na(refPop))     nR <- n1+n2
      if (refPop=="Cont"     & !is.na(refPop))     nR <- n3+n4
      if (refPop=="Opt"      & !is.na(refPop) & method.rr=="OP")   {
        p1.hat  <- sum(x1)/sum(n1)
        p2.hat  <- sum(x2)/sum(n2)
        p3.hat  <- sum(x3)/sum(n3)
        pa.hat  <- (p1.hat+p2.hat)/2
        pie.num <- (pa.hat*(1/n1 + 1/n2) + (0.01+(p1.hat-p2.hat)^2)/(n3*p3.hat))^(-1)
        pie     <- pie.num/sum(pie.num)
        nR      <- pie*sum(n1+n2+n3) }
      if (refPop=="Opt"      & !is.na(refPop) & method.rr=="TP")   {
        print("WARNING: The option refPop=Opt does not work for method.rr=TP; Defaults to refPop=All")
        nR <- n1+n2+n3+n4}


      ###Recalibrates x's and n's to reference group if needed
      if (!is.na(refPop)){
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
      if (length(goodStrata) < s & s > 1) print("WARNING: Dropping strata without any events (or events in treatment)")

      if (s > 1 & length(goodStrata) > 1) {x1 <- x1[goodStrata]
      x2 <- x2[goodStrata]
      x3 <- x3[goodStrata]
      x4 <- x4[goodStrata]
      n1 <- n1[goodStrata]
      n2 <- n2[goodStrata]
      n3 <- n3[goodStrata]
      n4 <- n4[goodStrata]}
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
      ncz    <- ncol(Z)
      badCol <- rep(0,ncz)
      for (i in 2:ncz) if (sum(Y[abs(Z[,i]) > 0])==0) badCol[i] <- 1
      Z <- Z[,badCol==0]
      init.b2 <- log(mean(Y[tr==0],na.rm=T))
      init.b  <- c(ifelse(is.na(init.b2) | init.b2 < -15, -15,init.b2),rep(0,ncol(Z)-1))
      init.k  <- mean(Y[tr!=0],na.rm=T)/exp(init.b[1])
      if (method.rr == "CO") init.k  <- ifelse(init.k-init.r/2 < 0.000001, init.r/2 + 0.000001, init.k )
      init.R2 <- mean(Y[tr==2],na.rm=T)/mean(Y[tr==0],na.rm=T)
      x1 <- x2 <- x3 <- x4 <- n1 <- n2 <- n3 <- n4 <- NA
    }



    EST <- LB <- UB <- NA
    if ( (sum(x3)>=2 & (sum(x4)>=2 | sum(n4)==0)) | substr(method.rr,1,2) == "CO" ){
      if (method.rr=="OP") EST <- getEst.OP(init.pa,init.t,init.p3,init.r,x1,x2,x3,n1,n2,n3,altParam,unCond)
      if (method.rr=="TP") EST <- getEst.TP(init.p4,init.tp,init.p3,init.r,x1,x2,x3,x4,n1,n2,n3,n4)
      if (method.rr=="CO") EST <- getEst.CO(init.k,init.b,init.r,Y,Z,tr)
      if (method.rr=="COR")EST <- getEst.COR(init.R2,init.b,init.R1,Y,Z,tr)
    }

    if ( (sum(x3)==0 | (sum(x4)==0 & sum(n4)>0)) &  substr(method.rr,1,2) != "CO")  EST     <- 0
    if (!is.na(EST) & LRT==0) LB    <- boundR(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,alpha2=alpha/2,type="LB",altParam=altParam,method.rr=method.rr,unCond=unCond,fast=fast)
    if (!is.na(EST) & LRT==0) UB    <- boundR(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,alpha2=alpha/2,type="UB",altParam=altParam,method.rr=method.rr,unCond=unCond,fast=fast)
    if (!is.na(EST) & LRT==1) try(LB <- boundR.LRT(init.r,x1,x2,x3,x4,n1,n2,n3,n4,Y,Z,tr,init.pa,init.p3,init.p4,init.t,init.tp,init.k,init.b,init.R1,init.R2,method.rr=method.rr,altParam=altParam,alpha2=alpha/2,type="LB"))
    if (!is.na(EST) & LRT==1) try(UB <- boundR.LRT(init.r,x1,x2,x3,x4,n1,n2,n3,n4,Y,Z,tr,init.pa,init.p3,init.p4,init.t,init.tp,init.k,init.b,init.R1,init.R2,method.rr=method.rr,altParam=altParam,alpha2=alpha/2,type="UB"))

    beta               <- NA
    if (substr(method.rr,1,2) == "CO"){
      cname              <- colnames(Z)
      if (length(cname)==0) cname <- rep("--",ncol(Z))
      bestEst            <- NA
      if (method.rr == "CO")  try(bestEst           <- getEst.CO( init.k,init.b,init.r,Y,Z,tr,allRes=1))
      if (method.rr == "COR") try(bestEst           <- getEst.COR(init.R2,init.b,init.R1,Y,Z,tr,allRes=1))
      if (length(bestEst) > 1){
        lbe               <- length(bestEst)
        if (method.rr == "CO")  III               <- fI.CO(bestEst[lbe], bestEst[1],bestEst[-c(1,lbe)],tr,length(Y),Z)
        if (method.rr == "COR") III               <- fI.COR(bestEst[lbe],bestEst[1],bestEst[-c(1,lbe)],tr,length(Y),Z)
        SE                <- rep(NA,lbe)
        try(SE            <- sqrt(diag(solve(III)))[-c(1,lbe)])
        betas             <- bestEst[-c(1,lbe)]
        beta              <- data.frame("name"=cname,"betas"=betas,"SE"=SE)
        colnames(beta)    <- c("","Beta","SE")
      }
    }

    out <- list("EST"=EST,"LB"=LB,"UB"=UB,"beta"=beta)

    ####end if deltaMethod==0
  }


  if (deltaMethod==1) out <- deltaMethodF(x1=x1,x2=x2,x3=x3,n1=n1,n2=n2,n3=n3,alpha=alpha,refPop=refPop)

  out

}


####################################################################################################################
#####
#####Delta Method
#####
####################################################################################################################



deltaMethodF <- function(x1=NA,x2=NA,x3=NA,n1=NA,n2=NA,n3=NA,alpha=0.05,refPop=NA){


  beta <- NA
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


  if (length(x1)>1 & (refPop != "Ind" | is.na(refPop)) ){

    if (is.na(refPop)) refPop <- "All"
    if (refPop=="All"      & !is.na(refPop))     nR <- n1+n2+n3
    if (refPop=="Treat"    & !is.na(refPop))     nR <- n1+n2
    if (refPop=="Cont"     & !is.na(refPop))     nR <- n3
    if (refPop=="Opt"      & !is.na(refPop))   {
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


  if (length(x1)>1 & refPop == "Ind"){


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

  list("EST"=EST,"LB"=LB,"UB"=UB,"beta"=beta)
}




####################################################################################################################
#####
#####Minimization routine runs nleqslv, but uses a wrapper to make the output look like rootSolve
#####
####################################################################################################################


nleqslv2 <- function(x,fn,...,positive=FALSE) {
  res <- nleqslv::nleqslv(x,fn,...)
  res$root   <- ifelse(positive==TRUE & res$x < 0, NA,res$x)
  res$f.root <- ifelse(positive==TRUE & res$x < 0, 10,res$fvec)
  res$estim.precis <- sum(abs(res$fvec))
  res}


####################################################################################################################
#####
#####Generates Data
#####
####################################################################################################################


###
### r = (p2 - p1)/p3
###
###Generates a Dataset (3 groups)
###input:
###r:                 difference in relative risks
###pa.s:              average risk in groups w/ drugs
###p3.s:              average risk in groups w/o drugs
###n1.s, n2.s, n3.s:  total number of subjects in each of the three groups (to be split among s strata)

genData.OP   <- function(pa.s,p3.s,r,n1.s,n2.s,n3.s,s=1){
  if (s>1)  tiltb <- seq(0.67,1.33,length.out=s)
  if (s==1) tiltb <- 1
  tilt  <- tiltb/sum(tiltb)
  r2   <- r/2
  pa   <- pa.s*tiltb
  p3   <- p3.s*tiltb
  p1   <- pa - r2*p3
  p2   <- pa + r2*p3
  n1   <- round(n1.s*tilt)
  n2   <- round(n2.s*rev(tilt))
  n3   <- round(rep(round(n3.s/s),s))
  x1   <- rbinom(s,n1,p1)
  x2   <- rbinom(s,n2,p2)
  x3   <- rbinom(s,n3,p3)
  list("x1"=x1,"x2"=x2,"x3"=x3,"n1"=n1,"n2"=n2,"n3"=n3)
}

###
### r = p2/p02 - p1/p01
### t = p2/p02 + p1/p01
###
###Generates a Dataset (4 groups; two control groups)
###input:
###r:                           difference in relative risks
###t:                           sum of vaccine efficacies
###p01.s, p02.s:                risk in each of the control groups
###n1.s, n2.s, n01.s, n02.s:    total number of subjects in each of the four groups (to be split among s strata)

genData.TP   <- function(t.s,r,p01.s,p02.s,n1.s,n2.s,n01.s,n02.s,s=1){
  t    <- rep(t.s,s)
  p1   <- rep(p01.s*(t-r)/2,s)
  p2   <- rep(p02.s*(t+r)/2,s)
  p01  <- rep(p01.s,s)
  p02  <- rep(p02.s,s)
  n1   <- rep(round(n1.s/s),s)
  n2   <- rep(round(n2.s/s),s)
  n01  <- rep(round(n01.s/s),s)
  n02  <- rep(round(n02.s/s),s)
  x1   <- rbinom(s,n1,p1)
  x2   <- rbinom(s,n2,p2)
  x01  <- rbinom(s,n01,p01)
  x02  <- rbinom(s,n02,p02)
  list("x1"=x1,"x2"=x2,"x01"=x01,"x02"=x02,"n1"=n1,"n2"=n2,"n01"=n01,"n02"=n02)
}

###
###p2 = exp(bZ)(k+r/2)
###p1 = exp(bZ)(k-r/2)
###p0 = exp(bZ)
###
###Generates the outcome (Y), the covariates (Z), the treatment status (tr) for the individuals
###input:
###n:         number of subjects
###r:         (see equation above)
###k:         (see equation above)
###b:         a set of coefficients for the covariates
###

genData.CO   <- function(r,k,b,n){
  ncov  <- length(b)-1
  npg   <- round(n/3)
  n2    <- 3*npg
  tr    <- rep(c(0,1,2),each=npg)
  Z     <- matrix(nrow=n2,ncol=ncov+1)
  Z[,1] <- matrix(rep(1,n2),ncol=1)
  for (i in 1:ncov) Z[,i+1]   <- rbinom(n2,1,0.5)
  p    <- exp(Z%*%b)*ifelse(tr==1,k-r/2,ifelse(tr==2,k+r/2,1))
  Y    <- rbinom(n2,1,p)
  list("Z"=Z,"Y"=Y,"tr"=tr)
}


####################################################################################################################
#####
#####Calculates the negative log-likelihood
#####
####################################################################################################################


###One Control Group:

fL.OP         <- function(t=NA,pa=NA,p3,r,x1,x2,x3,n1,n2,n3)   {
  s    <- length(p3)
  if (sum(!is.na(pa))>0) p1   <- pa   - r*p3/2
  if (sum(!is.na(pa))>0) p2   <- pa   + r*p3/2
  if (sum(!is.na(t))>0)  p1   <- t*p3 - r*p3/2
  if (sum(!is.na(t))>0)  p2   <- t*p3 + r*p3/2
  out  <- NA
  if (sum(p1>=0)==s & sum(p2>=0)==s & sum(p3>=0)==s &
      sum(p1<=1)==s & sum(p2<=1)==s & sum(p3<=1)==s) {p3 <- ifelse(p3<0.000001,0.000001,ifelse(p3==1,0.999999,p3))
      p1 <- ifelse(p1<0.000001,0.000001,ifelse(p1==1,0.999999,p1))
      p2 <- ifelse(p2<0.000001,0.000001,ifelse(p2==1,0.999999,p2))
      out2 <- -(x3*log(p3)+(n3-x3)*log(1-p3)+
                  x1*log(p1)+(n1-x1)*log(1-p1)+
                  x2*log(p2)+(n2-x2)*log(1-p2))
      out  <- sum(out2)}
  if (abs(out)>9999999999999 | is.na(out)) out <- 9999999999999
  out}



###Two Control groups

fL.TP         <- function(t,r,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02)   {
  s    <- length(p01)
  p1 <- p01*(t-r)/2
  p2 <- p02*(t+r)/2
  out  <- NA
  if (sum(p1>=0)==s & sum(p2>=0)==s & sum(p01>=0)==s & sum(p02>=0)==s &
      sum(p1<=1)==s & sum(p2<=1)==s & sum(p01<=1)==s & sum(p02<=1)==s) {p01 <- ifelse(p01<0.000001,0.000001,ifelse(p01==1,0.999999,p01))
      p02 <- ifelse(p02<0.000001,0.000001,ifelse(p02==1,0.999999,p02))
      p1  <- ifelse(p1<0.000001,0.000001,ifelse(p1==1,0.999999,p1))
      p2  <- ifelse(p2<0.000001,0.000001,ifelse(p2==1,0.999999,p2))
      out2 <- -(x02*log(p02)+(n02-x02)*log(1-p02)+
                  x01*log(p01)+(n01-x01)*log(1-p01)+
                  x1*log(p1)+(n1-x1)*log(1-p1)+
                  x2*log(p2)+(n2-x2)*log(1-p2))
      out  <- sum(out2)}
  if (abs(out)>9999999999999 | is.na(out)) out <- 9999999999999
  out}


###Individual-level Covarites

fL.CO         <- function(r,k,b,Y,Z,tr){
  p    <- exp(Z%*%b)*ifelse(tr==1,k-r/2,ifelse(tr==2,k+r/2,1))
  p    <- ifelse(p<=0.000001,0.000001,ifelse(p>=1,0.999999,p))
  out  <- -sum(Y*log(p)+(1-Y)*log(1-p))
  out
}


###Individual-level Covarites focused on relative risk

fL.COR        <- function(R1,R2,b,Y,Z,tr) fL.CO(r=(R2-R1),k=(R1+R2)/2,b,Y,Z,tr)




####################################################################################################################
#####
#####Calculates the score equations (i.e. derivative of the log-likelhood)
#####
####################################################################################################################



###One Control Group: Standard Parameterization

dL.OP1<- function(pa,p3,r,x1,x2,x3,n1,n2,n3)   {
  s   <- length(pa)
  p1=pa-r*p3/2
  p2=pa+r*p3/2
  p1=ifelse(x1>0 & p1<0.000001,0.000001,p1)
  p2=ifelse(x2>0 & p2<0.000001,0.000001,p2)
  p3=ifelse(x3>0 & p3<0.000001,0.000001,p3)
  p1d <- ifelse(p1==0,1,p1)
  p2d <- ifelse(p2==0,1,p2)
  p3d <- ifelse(p3==0,1,p3)
  out.r  <- sum(-x1*(p3/2)/(p1d) + (n1-x1)*(p3/2)/(1-p1) +
                  x2*(p3/2)/(p2d) - (n2-x2)*(p3/2)/(1-p2))
  out.pa <-      x1[1]/p1d[1] - (n1[1]-x1[1])/(1-p1[1])  +
    x2[1]/p2d[1] - (n2[1]-x2[1])/(1-p2[1])
  out.p3 <-      x3[1]/p3d[1] - (n3[1]-x3[1])/(1-p3[1])  -  x1[1]*(r/2)/(p1d[1]) + (n1[1]-x1[1])*(r/2)/(1-p1[1]) +
    x2[1]*(r/2)/(p2d[1]) - (n2[1]-x2[1])*(r/2)/(1-p2[1])
  out    <- c(out.r,out.pa,out.p3)
  if (s > 1) for (i in 2:s) {
    out.pa <-       x1[i]/p1d[i] - (n1[i]-x1[i])/(1-p1[i])  +
      x2[i]/p2d[i] - (n2[i]-x2[i])/(1-p2[i])
    out.p3 <-       x3[i]/p3d[i] - (n3[i]-x3[i])/(1-p3[i])  -  x1[i]*(r/2)/(p1d[i]) + (n1[i]-x1[i])*(r/2)/(1-p1[i]) +
      x2[i]*(r/2)/(p2d[i]) - (n2[i]-x2[i])*(r/2)/(1-p2[i])
    out    <- c(out,out.pa,out.p3)

  }
  out <- matrix(out,ncol=1)
  out
}

###One Control Group: Alternative Parameterization

dL.OP2<- function(t,p3,r,x1,x2,x3,n1,n2,n3)   {
  s   <- length(p3)
  p1=(t*p3-r*p3/2)
  p2=(t*p3+r*p3/2)
  p1=ifelse(x1>0 & p1<0.000001,0.000001,p1)
  p2=ifelse(x2>0 & p2<0.000001,0.000001,p2)
  p3=ifelse(x3>0 & p3<0.000001,0.000001,p3)
  p1d <- ifelse(p1==0,1,p1)
  p2d <- ifelse(p2==0,1,p2)
  p3d <- ifelse(p3==0,1,p3)
  out.r  <- sum(-x1*p3/(2*p1d) + (n1-x1)*p3/(2*(1-p1)) +
                  x2*p3/(2*p2d) - (n2-x2)*p3/(2*(1-p2)) )
  out.t <-  sum(x1*p3/p1d - (n1-x1)*p3/(1-p1) +
                  x2*p3/p2d - (n2-x2)*p3/(1-p2) )
  out.p3 <- x3[1]/p3d[1] - (n3[1]-x3[1])/(1-p3[1]) + x1[1]*(t-r/2)/p1d[1] - (n1[1]-x1[1])*(t-r/2)/(1-p1[1]) +
    x2[1]*(t+r/2)/p2d[1] - (n2[1]-x2[1])*(t+r/2)/(1-p2[1])
  out    <- c(out.r,out.t,out.p3)
  if (s > 1) for (i in 2:s) {
    out.p3 <-  x3[i]/p3d[i] - (n3[i]-x3[i])/(1-p3[i]) + x1[i]*(t-r/2)/p1d[i] - (n1[i]-x1[i])*(t-r/2)/(1-p1[i]) +
      x2[i]*(t+r/2)/p2d[i] - (n2[i]-x2[i])*(t+r/2)/(1-p2[i])
    out    <- c(out,out.p3)

  }
  out <- matrix(out,ncol=1)
  out
}


##Two Control Groups:

dL.TP    <- function(r,t,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02) {
  s <- length(p01)
  p1 <- p01*(t-r)/2
  p2 <- p02*(t+r)/2
  p01=ifelse(x01>0 & p01<0.000001,0.000001,p01)
  p02=ifelse(x02>0 & p02<0.000001,0.000001,p02)
  p1=ifelse(x1>0 & p1<0.000001,0.000001,p1)
  p2=ifelse(x2>0 & p2<0.000001,0.000001,p2)
  q1 <- 1-p1
  q2 <- 1-p2
  q01<- 1-p01
  q02<- 1-p02
  out.r   <- sum(-x1/(t-r)+(n1-x1)*p01/(2*q1)+x2/(t+r)-(n2-x2)*p02/(2*q2))
  out.t   <- x1[1]/(t[1]-r)-(n1[1]-x1[1])*p01[1]/(2*q1[1])+x2[1]/(t[1]+r)-(n2[1]-x2[1])*p02[1]/(2*q2[1])
  out.P01 <- x01[1]/p01[1]-(n01[1]-x01[1])/q01[1]+x1[1]/p01[1]-(n1[1]-x1[1])*(t[1]-r)/(2*q1[1])
  out.P02 <- x02[1]/p02[1]-(n02[1]-x02[1])/q02[1]+x2[1]/p02[1]-(n2[1]-x2[1])*(t[1]+r)/(2*q2[1])
  out <- c(out.r,out.t,out.P01,out.P02)
  if (s > 1) for (i in 2:s){
    out     <- c(out, x1[i]/(t[i]-r)-(n1[i]-x1[i])*p01[i]/(2*q1[i])+x2[i]/(t[i]+r)-(n2[i]-x2[i])*p02[i]/(2*q2[i]),
                 x01[i]/p01[i]-(n01[i]-x01[i])/q01[i]+x1[i]/p01[i]-(n1[i]-x1[i])*(t[i]-r)/(2*q1[i]),
                 x02[i]/p02[i]-(n02[i]-x02[i])/q02[i]+x2[i]/p02[i]-(n2[i]-x2[i])*(t[i]+r)/(2*q2[i]))
  }
  matrix(out,ncol=1)
}


####Individual-Level covariates

dL.CO    <- function(r,k,b,tr,n,Z,Y) {
  nc     <- length(b)
  ezb    <- exp(Z%*%b)
  p      <- ezb*ifelse(tr==1,k-r/2,ifelse(tr==2,k+r/2,1))
  p      <- ifelse(p<0.000001,0.000001,p)
  q      <- 1-p
  q      <- ifelse(q<0.000001,0.000001,q)
  out.r  <- sum((-Y/(2*(k-r/2)) + (1-Y)*ezb /(2*q))*ifelse(tr==1,1,0) +   (Y/(2*(k+r/2)) - (1-Y)*ezb/(2*q))*ifelse(tr==2,1,0))
  out.k  <- sum( (Y/(  (k-r/2)) - (1-Y)*ezb /(  q))*ifelse(tr==1,1,0) +   (Y/(  (k+r/2)) - (1-Y)*ezb/(  q))*ifelse(tr==2,1,0))
  out    <- c(out.r,out.k)
  for (i in 1:nc) out <- c(out,sum(Y*Z[,i]-(1-Y)*Z[,i]*p/q))
  matrix(out,ncol=1)
}



####Individual-Level covariates: relative risk

dL.COR    <- function(R1,R2,b,tr,n,Z,Y) {
  nc     <- length(b)
  ezb    <- exp(Z%*%b)
  p      <- ezb*ifelse(tr==1,R1,ifelse(tr==2,R2,1))
  p      <- ifelse(p<0.000001,0.000001,p)
  q      <- 1-p
  q      <- ifelse(q<0.000001,0.000001,q)
  out.R1  <- sum( (Y/R1 - (1-Y)*ezb /q)*ifelse(tr==1,1,0))
  out.R2  <- sum( (Y/R2 - (1-Y)*ezb /q)*ifelse(tr==2,1,0))
  out    <- c(out.R1,out.R2)
  for (i in 1:nc) out <- c(out,sum(Y*Z[,i]-(1-Y)*Z[,i]*p/q))
  matrix(out,ncol=1)
}



####################################################################################################################
#####
#####Calculates the negative of the information matrices
#####
####################################################################################################################


###One Control Group: Standard Parameterization

fI.OP1   <- function(pa,p3,r,n1,n2,n3){

  r2 <- r/2
  p1 <- pa - r2*p3
  p2 <- pa + r2*p3
  if (sum(p1==0) > 0 | sum(p2==0) > 0) { pa <- ifelse(p1==0 | p2==0,pa+0.000001,pa)
  p1 <- pa - r2*p3
  p2 <- pa + r2*p3
  }
  q1 <- 1-p1
  q2 <- 1-p2
  q3 <- 1-p3
  #p1      <- ifelse(p1>=0 & p1 < 0.000001 ,0.000001,p1)
  #q1      <- ifelse(q1>=0 & q1 < 0.000001, 0.000001,q1)
  #p2      <- ifelse(p2>=0 & p2 < 0.000001 ,0.000001,p2)
  #q2      <- ifelse(q2>=0 & q2 < 0.000001, 0.000001,q2)
  #p3      <- ifelse(p3>=0 & p3 < 0.000001 ,0.000001,p3)
  #q3      <- ifelse(q3>=0 & q3 < 0.000001, 0.000001,q3)
  s  <- length(pa)
  prob <- 0
  if(sum(p1<0) > 0 | sum(p1>1) > 0 | sum(p2<0) > 0 | sum(p2>1) > 0) prob <- 1
  sv <- ifelse(prob==0,0,NA)
  I      <- matrix(sv,nrow=(2*s+1),ncol=(2*s+1))
  if (prob==0){
    I[1,1] <- dL2.RR   <- sum(-(n1*(p3/2)^2)/(p1*q1) - (n2*(p3/2)^2)/(p2*q2))
    for (i in 1:s){
      sp <- 2+2*(i-1)
      ep <- sp + 1
      I[ep,ep] <- dL2.P3P3 <- -n3[i]/(p3[i]*q3[i]) - n1[i]*r2^2/(p1[i]*q1[i]) - n2[i]*r2^2/(p2[i]*q2[i])
      I[sp,sp] <- dL2.PAPA <- -n1[i]/(p1[i]*q1[i]) - n2[i]/(p2[i]*q2[i])

      I[sp,ep] <- I[ep,sp] <- dL2.P3PA <-  n1[i]*r2/(p1[i]*q1[i]) - n2[i]*r2/(p2[i]*q2[i])
      I[1,sp]  <- I[sp,1]  <- dL2.PAR  <-  n1[i]*(p3[i]/2)/(p1[i]*q1[i]) - n2[i]*(p3[i]/2)/(p2[i]*q2[i])
      I[1,ep]  <- I[ep,1]  <- dL2.P3R  <- -n1[i]*(pa[i]/2)/p1[i] +
        n1[i]*((1-pa[i])/2)/q1[i] +
        n2[i]*(pa[i]/2)/p2[i] +
        -n2[i]*((1-pa[i])/2)/q2[i]
    }
  }
  return(-I)
}


###One Control Group: Alternative Parameterization

fI.OP2   <- function(t,p3,r,n1,n2,n3){

  r2 <- r/2
  p1 <- t*p3 - r2*p3
  p2 <- t*p3 + r2*p3
  q1 <- 1-p1
  q2 <- 1-p2
  q3 <- 1-p3
  p1      <- ifelse(p1>=0 & p1 < 0.000001 ,0.000001,p1)
  q1      <- ifelse(q1>=0 & q1 < 0.000001, 0.000001,q1)
  p2      <- ifelse(p2>=0 & p2 < 0.000001 ,0.000001,p2)
  q2      <- ifelse(q2>=0 & q2 < 0.000001, 0.000001,q2)
  p3      <- ifelse(p3>=0 & p3 < 0.000001 ,0.000001,p3)
  q3      <- ifelse(q3>=0 & q3 < 0.000001, 0.000001,q3)
  s  <- length(p3)
  prob <- 0
  if(sum(p1<0) > 0 | sum(p1>1) > 0 | sum(p2<0) > 0 | sum(p2>1) > 0) prob <- 1
  sv <- ifelse(prob==0,0,NA)
  I      <- matrix(sv,nrow=(s+2),ncol=(s+2))
  if (prob==0){
    I[1,1] <- dL2.RR   <- sum(-(n1*(p3/2)^2)/(p1*q1) - (n2*(p3/2)^2)/(p2*q2))
    I[2,2] <- dL2.TT   <- sum(-(n1*p3^2)/(p1*q1) - (n2*p3^2)/(p2*q2))
    I[1,2] <- I[2,1]   <- dL2.TR  <-  sum(n1*p3^2/(2*p1*q1) - n2*p3^2/(2*p2*q2))

    for (i in 1:s){
      sp <- 2+i
      I[sp,sp] <- dL2.P3P3 <- -n3[i]/(p3[i]*q3[i]) - n1[i]*(t-r2)^2/(p1[i]*q1[i]) - n2[i]*(t+r2)^2/(p2[i]*q2[i])
      I[1,sp]  <- I[sp,1]  <- dL2.P3R  <-   n1[i]/(2*q1[i]) - n2[i]/(2*q2[i])
      I[2,sp]  <- I[sp,2]  <- dL2.P3T  <-  -n1[i]/q1[i]     - n2[i]/q2[i]
    }
  }
  return(-I)
}

###Two Control Groups:

fI.TP <- function(t,r,p01,p02,n1,n2,n01,n02) {
  p1 <- p01*(t-r)/2
  p2 <- p02*(t+r)/2
  q1 <- 1-p1
  q2 <- 1-p2
  q01<- 1-p01
  q02<- 1-p02
  s  <- length(p1)
  prob <- 0
  if(sum(p1<0) > 0 | sum(p1>1) > 0 | sum(p2<0) > 0 | sum(p2>1) > 0 | sum(p01<0) > 0 | sum(p01>1) > 0 | sum(p02<0) > 0 | sum(p02>1) > 0) prob <- 1
  sv <- ifelse(prob==0,0,NA)
  I      <- matrix(sv,nrow=(3*s+1),ncol=(3*s+1))
  I[1,1]     <- sum(-n1*p01^2/(4*p1*q1) -n2*p02^2/(4*p2*q2))
  I[2,2]     <-  -n1[1]*p01[1]^2/(4*p1[1]*q1[1]) -n2[1]*p02[1]^2/(4*p2[1]*q2[1])
  I[3,3]     <- out.P01P01 <- -n01[1]/(p01[1]*q01[1]) - (n1[1]*(t[1]-r)^2)/(4*p1[1]*q1[1])
  I[4,4]     <- out.P02P02 <- -n02[1]/(p02[1]*q02[1]) - (n1[1]*(t[1]+r)^2)/(4*p2[1]*q2[1])
  I[1,3]     <- I[3,1]     <- out.P01r  <- n1[1]/(2*q1[1])
  I[1,4]     <- I[1,4]     <- out.P02r <- -n2[1]/(2*q2[1])
  I[2,3]     <- I[3,2]     <- out.P01t <- -n1[1]/(2*q1[1])
  I[2,4]     <- I[4,2]     <- out.P02t <- -n2[1]/(2*q2[1])
  I[1,2]     <- I[2,1]     <- out.tr     <- n1[1]*p01[1]^2/(4*p1[1]*q1[1]) -n2[1]*p02[1]^2/(4*p2[1]*q2[1])
  if (s > 1) for (i in 2:s) {
    sv <- (i-1)*3+2
    I[sv,sv]     <-  -n1[i]*p01[i]^2/(4*p1[i]*q1[i]) -n2[i]*p02[i]^2/(4*p2[i]*q2[i])
    I[sv+1,sv+1]     <- out.P01P01 <- -n01[i]/(p01[i]*q01[i]) - (n1[i]*(t[i]-r)^2)/(4*p1[i]*q1[i])
    I[sv+2,sv+2]     <- out.P02P02 <- -n02[i]/(p02[i]*q02[i]) - (n1[i]*(t[i]+r)^2)/(4*p2[i]*q2[i])
    I[1,sv+1]     <- I[sv+1,1]     <- out.P01r  <- n1[i]/(2*q1[i])
    I[1,sv+2]     <- I[sv+2.1]     <- out.P02r <- -n2[i]/(2*q2[i])
    I[2,sv+1]     <- I[sv+1,2]     <- out.P01t <- -n1[i]/(2*q1[i])
    I[2,sv+2]     <- I[sv+2,2]     <- out.P02t <- -n2[i]/(2*q2[i])
    I[1,sv]      <- I[sv,1]     <- out.tr     <- n1[i]*p01[i]^2/(4*p1[i]*q1[i]) -n2[i]*p02[i]^2/(4*p2[i]*q2[i])

  }
  -I
}


###Individual-level covariates

fI.CO <- function(r,k,b,tr,n,Z) {
  nc     <- length(b)
  ezb    <- exp(Z%*%b)
  #ezb    <- ifelse(tr==1 & ezb*(k-r/2) < 0.000001,0.000001/(k-r/2),
  #                 ifelse(tr==2 & ezb*(k+r/2) < 0.000001,0.000001/(k+r/2),ezb))
  ezb2   <- ezb^2
  p      <- ezb*ifelse(tr==1,k-r/2,ifelse(tr==2,k+r/2,1))
  q      <- 1-p
  p      <- ifelse(p>=0 & p < 0.000001 ,0.000001,p)
  q      <- ifelse(q>=0 & q < 0.000001, 0.000001,q)
  Zmat   <- matrix(nrow=nc,ncol=nc)
  I      <- matrix(nrow=(nc+2),ncol=(nc+2))
  for (i in 1:nc) for (j in 1:nc)  Zmat[i,j] <- sum(-Z[,i]*Z[,j]*p/q)
  I[-c(1,2),-c(1:2)] <- Zmat
  I[2,2]  <-           sum((-p/(  (k-r/2)^2) - ezb2/   q )*ifelse(tr==1,1,0) + (-p/(  (k+r/2)^2) - ezb2/   q )*ifelse(tr==2,1,0))
  I[1,1]  <-           sum((-p/(4*(k-r/2)^2) - ezb2/(4*q))*ifelse(tr==1,1,0) + (-p/(4*(k+r/2)^2) - ezb2/(4*q))*ifelse(tr==2,1,0))
  I[1,2]  <- I[2,1] <- sum(( p/(2*(k-r/2)^2) + ezb2/(2*q))*ifelse(tr==1,1,0) + (-p/(2*(k+r/2)^2) - ezb2/(2*q))*ifelse(tr==2,1,0))
  for (i in 1:nc) I[1,i+2] <- I[i+2,1] <- sum(  (Z[,i]*ezb/(2*q))*ifelse(tr==1,1,0) - (Z[,i]*ezb/(2*q))*ifelse(tr==2,1,0))
  for (i in 1:nc) I[2,i+2] <- I[i+2,2] <- -sum( (Z[,i]*ezb/(  q))*ifelse(tr==1,1,0) + (Z[,i]*ezb/(q))*ifelse(tr==2,1,0))
  -I
}






###Individual-level covariates: Relative Risk


fI.COR <- function(R1,R2,b,tr,n,Z) {
  nc     <- length(b)
  ezb    <- exp(Z%*%b)
  ezb2   <- ezb^2
  p      <- ezb*ifelse(tr==1,R1,ifelse(tr==2,R2,1))
  q      <- 1-p
  p      <- ifelse(p>=0 & p < 0.000001 ,0.000001,p)
  q      <- ifelse(q>=0 & q < 0.000001, 0.000001,q)
  Zmat   <- matrix(nrow=nc,ncol=nc)
  I      <- matrix(nrow=(nc+2),ncol=(nc+2))
  for (i in 1:nc) for (j in 1:nc)  Zmat[i,j] <- sum(-Z[,i]*Z[,j]*p/q)
  I[-c(1,2),-c(1:2)] <- Zmat
  I[2,2]  <- sum((-p/(R2)^2 - ezb2/q)*ifelse(tr==2,1,0))
  I[1,1]  <- sum((-p/(R1)^2 - ezb2/q)*ifelse(tr==1,1,0))
  I[1,2]  <- I[2,1] <- 0
  for (i in 1:nc) I[1,i+2] <- I[i+2,1] <- -sum(  (Z[,i]*ezb/(q))*ifelse(tr==1,1,0))
  for (i in 1:nc) I[2,i+2] <- I[i+2,2] <- -sum(  (Z[,i]*ezb/(q))*ifelse(tr==2,1,0))
  -I
}


####################################################################################################################
#####
#####Funtions that calculate the log-likelihood or score vectors ... using different formats of input
#####
####################################################################################################################


fL.pap3r         <- function(pap3r,x1,x2,x3,n1,n2,n3,deriv=0) {s <- length(x1)
pa<- pap3r[1:s]
p3<- pap3r[(s+1):(2*s)]
r <- pap3r[2*s+1]
if (deriv==0) out <- fL.OP(NA,pa,p3,r,x1,x2,x3,n1,n2,n3)
if (deriv==1) out <- dL.OP1(pa,p3,r,x1,x2,x3,n1,n2,n3)
out
}


fL.pap3         <- function(pap3,x1,x2,x3,n1,n2,n3,r,deriv=0) {s <- length(x1)
pa<- pap3[1:s]
p3<- pap3[(s+1):(2*s)]
r2=r
if (deriv==0) out <- fL.OP(NA,pa,p3,r2,x1,x2,x3,n1,n2,n3)
if (deriv==1) out <- dL.OP1(pa,p3,r2,x1,x2,x3,n1,n2,n3)[-1]
out
}


fL.pap3.parms <- function(pap3,parms) { s     <- parms[1]
x1    <- parms[2:(s+1)]
x2    <- parms[(s+2):(2*s+1)]
x3    <- parms[(2*s+2):(3*s+1)]
n1    <- parms[(3*s+2):(4*s+1)]
n2    <- parms[(4*s+2):(5*s+1)]
n3    <- parms[(5*s+2):(6*s+1)]
r     <- parms[(6*s+2)]
deriv <- parms[(6*s+3)]
fL.pap3(pap3,x1=x1,x2=x2,x3=x3,n1=n1,n2=n2,n3=n3,r=r,deriv=deriv)
}


fL.pap3r.parms <- function(pap3r,parms) { s     <- parms[1]
x1    <- parms[2:(s+1)]
x2    <- parms[(s+2):(2*s+1)]
x3    <- parms[(2*s+2):(3*s+1)]
n1    <- parms[(3*s+2):(4*s+1)]
n2    <- parms[(4*s+2):(5*s+1)]
n3    <- parms[(5*s+2):(6*s+1)]
deriv <- parms[(6*s+2)]
fL.pap3r(pap3r,x1=x1,x2=x2,x3=x3,n1=n1,n2=n2,n3=n3,deriv=deriv)
}

fL.pap3r.parms.m5 <- function(pap3r,parms) fL.pap3r.parms(pap3r-c(rep(0,2*parms[1]),5) ,parms)


fL.tp3r         <- function(tp3r,x1,x2,x3,n1,n2,n3,deriv=0) {s <- length(x1)
t <- tp3r[1]
p3<- tp3r[2:(s+1)]
r <- tp3r[s+2]
if (deriv==0) out=fL.OP(t,NA,p3,r,x1,x2,x3,n1,n2,n3)
if (deriv==1) out=dL.OP2(t,p3,r,x1,x2,x3,n1,n2,n3)
out
}

fL.tp3         <- function(tp3,x1,x2,x3,n1,n2,n3,r,deriv=0) {s <- length(x1)
t <- tp3[1]
p3<- tp3[2:(s+1)]
if (deriv==0) out=fL.OP(t,NA,p3,r,x1,x2,x3,n1,n2,n3)
if (deriv==1) out=dL.OP2(t,p3,r,x1,x2,x3,n1,n2,n3)[-1]
out
}


fL.tp3.parms <- function(tp3,parms) { s     <- parms[1]
x1    <- parms[2:(s+1)]
x2    <- parms[(s+2):(2*s+1)]
x3    <- parms[(2*s+2):(3*s+1)]
n1    <- parms[(3*s+2):(4*s+1)]
n2    <- parms[(4*s+2):(5*s+1)]
n3    <- parms[(5*s+2):(6*s+1)]
r     <- parms[(6*s+2)]
deriv <- parms[(6*s+3)]
fL.tp3(tp3,x1=x1,x2=x2,x3=x3,n1=n1,n2=n2,n3=n3,r=r,deriv=deriv)
}


fL.tp3r.parms <- function(tp3r,parms) { s     <- parms[1]
x1    <- parms[2:(s+1)]
x2    <- parms[(s+2):(2*s+1)]
x3    <- parms[(2*s+2):(3*s+1)]
n1    <- parms[(3*s+2):(4*s+1)]
n2    <- parms[(4*s+2):(5*s+1)]
n3    <- parms[(5*s+2):(6*s+1)]
deriv <- parms[(6*s+2)]
fL.tp3r(tp3r,x1=x1,x2=x2,x3=x3,n1=n1,n2=n2,n3=n3,deriv=deriv)
}
fL.tp3r.parms.m5 <- function(tp3r,parms) fL.tp3r.parms(tp3r-c(rep(0,1+parms[1]),5) ,parms)



fL.tp01p02r         <- function(tp01p02r,x1,x2,x01,x02,n1,n2,n01,n02,deriv=0) {s <- length(x1)
t  <- tp01p02r[1:s]
p01 <- tp01p02r[(s+1):(2*s)]
p02 <- tp01p02r[(2*s+1):(3*s)]
r <- tp01p02r[3*s+1]
if (deriv==0) out <- fL.TP(t,r,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02)
if (deriv==1) out <- dL.TP(r,t,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02)
out
}


fL.tp01p02         <- function(tp01p02,x1,x2,x01,x02,n1,n2,n01,n02,r,deriv=0) {s <- length(x1)
t  <- tp01p02[1:s]
p01 <- tp01p02[(s+1):(2*s)]
p02 <- tp01p02[(2*s+1):(3*s)]
if (deriv==0) out <- fL.TP(t,r,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02)
if (deriv==1) out <- dL.TP(r,t,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02)[-1]
out
}


fL.tp01p02.parms <- function(tp01p02,parms) { s     <- parms[1]
x1    <- parms[2:(s+1)]
x2    <- parms[(s+2):(2*s+1)]
x01    <- parms[(2*s+2):(3*s+1)]
x02    <- parms[(3*s+2):(4*s+1)]
n1    <- parms[(4*s+2):(5*s+1)]
n2    <- parms[(5*s+2):(6*s+1)]
n01    <- parms[(6*s+2):(7*s+1)]
n02    <- parms[(7*s+2):(8*s+1)]
r      <- parms[(8*s+2)]
deriv <- parms[(8*s+3)]
fL.tp01p02(tp01p02,x1,x2,x01,x02,n1,n2,n01,n02,r,deriv)
}


fL.tp01p02r.parms <- function(tp01p02r,parms) {  s     <- parms[1]
x1    <- parms[2:(s+1)]
x2    <- parms[(s+2):(2*s+1)]
x01    <- parms[(2*s+2):(3*s+1)]
x02    <- parms[(3*s+2):(4*s+1)]
n1    <- parms[(4*s+2):(5*s+1)]
n2    <- parms[(5*s+2):(6*s+1)]
n01    <- parms[(6*s+2):(7*s+1)]
n02    <- parms[(7*s+2):(8*s+1)]
deriv <- parms[(8*s+2)]
fL.tp01p02r(tp01p02r,x1,x2,x01,x02,n1,n2,n01,n02,deriv)
}
fL.tp01p02r.parms.m5 <- function(tp01p02r,parms) fL.tp01p02r.parms(tp01p02r-c(rep(0,3*parms[1]),5) ,parms)



fL.kbr         <- function(kbr,Y,Z,tr,deriv=0){
  k <- kbr[1]
  b <- kbr[2:(1+ncol(Z))]
  r <- kbr[(2+ncol(Z))]
  if (deriv==0) out <- fL.CO(r,k,b,Y,Z,tr)
  if (deriv==1) out <- dL.CO(r,k,b,tr,n=length(Y),Z,Y)
  out
}

fL.kb         <- function(kb,Y,Z,tr,r,deriv=0){
  k <- kb[1]
  b <- kb[2:(1+ncol(Z))]
  if (deriv==0) out <- fL.CO(r,k,b,Y,Z,tr)
  if (deriv==1) out <- dL.CO(r,k,b,tr,n=length(Y),Z,Y)[-1]
  out
}


fL.kb.parms <- function(kb,parms) {   n     <- parms[1]
r     <- parms[2]
deriv <- parms[3]
Y     <- parms[(4:(n+3))]
tr    <- parms[(n+4):(2*n+3)]
Z     <- matrix(parms[-c(1:(2*n+3))],byrow=F,nrow=n)
fL.kb(kb,Y,Z,tr,r,deriv)
}


fL.kbr.parms <- function(kbr,parms) {   n     <- parms[1]
deriv <- parms[2]
Y     <- parms[(3:(n+2))]
tr    <- parms[(n+3):(2*n+2)]
Z     <- matrix(parms[-c(1:(2*n+2))],byrow=F,nrow=n)
fL.kbr(kbr,Y,Z,tr,deriv)
}



#########################################################




fL.kbr2         <- function(kbr,Y,Z,tr,deriv=0){
  k <- kbr[1]
  b <- kbr[2:(1+ncol(Z))]
  r <- kbr[(2+ncol(Z))]
  if (deriv==0) out <- fL.COR(r,k,b,Y,Z,tr)
  if (deriv==1) out <- dL.COR(r,k,b,tr,n=length(Y),Z,Y)
  out
}

fL.kb2         <- function(kb,Y,Z,tr,r,deriv=0){
  k <- kb[1]
  b <- kb[2:(1+ncol(Z))]
  if (deriv==0) out <- fL.COR(r,k,b,Y,Z,tr)
  if (deriv==1) out <- dL.COR(r,k,b,tr,n=length(Y),Z,Y)[-1]
  out
}


fL.kb.parms2 <- function(kb,parms) {   n     <- parms[1]
r     <- parms[2]
deriv <- parms[3]
Y     <- parms[(4:(n+3))]
tr    <- parms[(n+4):(2*n+3)]
Z     <- matrix(parms[-c(1:(2*n+3))],byrow=F,nrow=n)
fL.kb2(kb,Y,Z,tr,r,deriv)
}


fL.kbr.parms2 <- function(kbr,parms) {   n     <- parms[1]
deriv <- parms[2]
Y     <- parms[(3:(n+2))]
tr    <- parms[(n+3):(2*n+2)]
Z     <- matrix(parms[-c(1:(2*n+2))],byrow=F,nrow=n)
fL.kbr2(kbr,Y,Z,tr,deriv)
}





####################################################################################################################
#####
##### Identifies teh conditional and unconditional MLE
#####
####################################################################################################################

optim.range <- function(pa=NA,p3=NA,t=NA,r=NA,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="OP",x1,x2,x3,x4,n1,n2,n3,n4,n,Z,Y,tr,fix.r=0,allRes=0,s=NA,allParam=NA,paramType=NA,nparam=NA,curIndex=NA,curParam=NA){

  if (method.rr=="OP" & altParam==0)  strata.s        <- c(rep(c(1:s),2),0)
  if (method.rr=="OP" & altParam==1)  strata.s        <- c(0,c(1:s),0)


  if (method.rr=="OP" & altParam==0){ if (curParam == "pa")  {min.v <-   max(  abs(allParam[paramType=="r"]) * allParam[paramType=="p3" & strata.s==strata.s[curIndex]]/2,0.000001)
  max.v <-   min(1-abs(allParam[paramType=="r"]) * allParam[paramType=="p3" & strata.s==strata.s[curIndex]]/2,0.999999)}
    if (curParam == "r")                                     {min.v <-   max(c(-2*(1-allParam[paramType =="pa"])/allParam[paramType=="p3"],-2*allParam[paramType =="pa"]/allParam[paramType =="p3"]))
    max.v <-   min(c(2*allParam[paramType =="pa"]/allParam[paramType=="p3"],2*(1-allParam[paramType =="pa"])/allParam[paramType =="p3"])) }
    if (curParam == "p3")  {paJ <- allParam[paramType =="pa" & strata.s==strata.s[curIndex]]
    rJ  <- allParam[paramType=="r"]
    if (rJ != 0) max.v <-   min(c(2*(1-paJ)/abs(rJ),2*paJ/abs(rJ),  0.99999))
    min.v <-  0.00001
    if (rJ == 0) max.v <-  0.99999
    }}



  if (method.rr=="OP" & altParam==1){ if (curParam == "t")  { min.v <-   max(      c(abs(allParam[paramType=="r"])/2) )
  max.v <-   min(      c((1 + allParam[paramType=="r"] * allParam[paramType=="p3"]/2)/allParam[paramType=="p3"],
                         (1 - allParam[paramType=="r"] * allParam[paramType=="p3"]/2)/allParam[paramType=="p3"]))}
    if (curParam == "r")   {min.v <-   max( c( -2*(1-allParam[paramType =="t"]*allParam[paramType=="p3"])/allParam[paramType=="p3"], -2*allParam[paramType =="t"]    ))
    max.v <-   min( c( 2*(1-allParam[paramType =="t"]*allParam[paramType=="p3"])/allParam[paramType=="p3"], 2*allParam[paramType =="t"]    )) }
    if (curParam == "p3")  {min.v <-   0.000001
    max.v <-   min( 1/(allParam[paramType =="t"]+abs(allParam[paramType=="r"])/2),0.99999)}}


  if (method.rr=="CO")              { if (curParam == "k")  {  rrr   <- allParam[length(allParam)]
  min.v <- 1.0001*abs(rrr/2) + 0.000001*ifelse(abs(rrr) < 0.000001,1,0)
  max.v <- min(1/max(exp(Z%*%b))-rrr/2,1/max(exp(Z%*%b))+rrr/2,3)}
    if (curParam == "b")  {  rrr  <- allParam[length(allParam)]
    kkk  <- allParam[1]
    bbb  <- allParam[-c(1,curIndex,length(allParam))]
    ppp  <- exp(Z[,-(curIndex-1)]%*%matrix(bbb,ncol=1))*ifelse(tr==1,kkk-rrr/2,ifelse(tr==2,kkk+rrr/2,1))
    ppp  <- ifelse(ppp<0.000001,0.000001,ppp)
    min.v<- max(  ifelse(Z[,(curIndex-1)]<0,log(1/ppp)/Z[,(curIndex-1)],-10^10) )
    max.v<- min(  ifelse(Z[,(curIndex-1)]>0,log(1/ppp)/Z[,(curIndex-1)], 10^10) )}
    if (curParam == "r")  {  kkk  <- allParam[1]
    bbb  <- allParam[-c(1,length(allParam))]
    min.v <- min(ifelse(tr==2,-2*kkk,
                        ifelse(tr==1,-2*(1/exp(Z%*%matrix(bbb,ncol=1))-kkk),-5)))
    max.v <- min(ifelse(tr==1,2*kkk,
                        ifelse(tr==1, 2*(1/exp(Z%*%matrix(bbb,ncol=1))-kkk),5)))}}



  if (method.rr=="COR")             {if (curParam == "R2" | curParam == "R1") {min.v <- 0.000001
  max.v <- 1/allParam[paramType=="p3"]
  mult2 <- 0.1}
    if (curParam == "b")                     {R2R2  <- allParam[length(allParam)]
    R1R1  <- allParam[1]
    bbb   <- allParam[-c(1,curIndex,length(allParam))]
    ppp   <- exp(Z[,-(curIndex-1)]%*%matrix(bbb,ncol=1))*ifelse(tr==1,R1R1,ifelse(tr==2,R2R2,1))
    min.v <- max(  ifelse(Z[,(curIndex-1)]<0,log(1/ppp)/Z[,(curIndex-1)],-10^10) )
    max.v  <- min(  ifelse(Z[,(curIndex-1)]>0,log(1/ppp)/Z[,(curIndex-1)],10^10) )}
    mult2  <- 2}


  list("min.v"=min.v,"max.v"=max.v)
}


optim.force <-  function(pa=NA,p3=NA,t=NA,r=NA,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="OP",x1,x2,x3,x4,n1,n2,n3,n4,n,Z,Y,tr,fix.r=0,allRes=0){
  s              <- length(x1)
  if (method.rr=="OP" & altParam==0)  allParam       <- c(pa,p3,r)
  if (method.rr=="OP" & altParam==1)  allParam       <- c(t,p3,r)
  if (method.rr=="CO")                allParam       <- c(k,b,r)
  if (method.rr=="COR")               allParam       <- c(R2,b,R1)
  if (method.rr=="TP")                allParam       <- c(tp,p3,p4,r)


  paramType <- rep(NA,length(allParam))
  if (method.rr=="OP" & altParam==0)  paramType       <- c(rep("pa",length(pa)),rep("p3",length(p3)),rep("r",length(r)))

  if (method.rr=="OP" & altParam==1)  paramType       <- c(rep("t",length(t)),  rep("p3",length(p3)),rep("r",length(r)))

  if (method.rr=="CO")                paramType       <- c(rep("k",length(k)),  rep("b",length(b)),  rep("r",length(r)))
  if (method.rr=="COR")               paramType       <- c(rep("R2",length(R2)),rep("b",length(b)),  rep("R1",length(R1)))
  if (method.rr=="TP")                paramType       <- c(rep("tp",length(pa)),rep("p3",length(p3)),rep("p4",length(p4)),rep("r",length(r)))


  if (fix.r==0)  nparam <- length(allParam)
  if (fix.r==1)  nparam <- length(allParam)-1

  prevMin <- 99999999999
  mult    <- 1
  out     <- NA
  #print("WARNING: The main function is now identifying the MLE for the liklihood by searching over a grid. This is a time-intensive step.")
  #print("WARNING: If there are a large number of strata (or covariates), this grid search may take > 10-20 minutes.")
  for (nattempt in 1:10000){
    for (curIndex in 1:nparam){

      curParam <- paramType[curIndex]
      range  <- optim.range(pa=pa,p3=p3,t=t,r=r,k=k,b=b,p4=p4,tp=tp,R1=R1,R2=R2,altParam=altParam,method.rr=method.rr,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,n=n,Z=Z,Y=Y,tr=tr,fix.r=fix.r,allRes=allRes,s=s,allParam=allParam,paramType=paramType,nparam=nparam,curIndex=curIndex,curParam=curParam)
      min.v  <- range$min.v
      max.v  <- range$max.v


      if (method.rr=="OP") if(allParam[curIndex] >= min.v & allParam[curIndex] <= max.v)  poss.vals <- sort(c(allParam[curIndex],seq(max(allParam[curIndex]-0.1*mult,min.v),min(allParam[curIndex]+0.1*mult,max.v),length.out=10)))
      if (method.rr=="OP") if((allParam[curIndex] < min.v | allParam[curIndex] > max.v))  poss.vals <- sort(c(seq(max(allParam[curIndex]-0.1*mult,min.v),min(allParam[curIndex]+0.1*mult,max.v),length.out=10)))



      if (method.rr=="CO") if( allParam[curIndex]  >= min.v & allParam[curIndex] <= max.v) poss.vals <- sort(c(   allParam[curIndex]      ,seq(max(allParam[curIndex]-2*mult,  min.v),min(allParam[curIndex]+2*mult,  max.v),length.out=10)))
      if (method.rr=="CO") if( (allParam[curIndex] < min.v  | allParam[curIndex] > max.v)) poss.vals <- sort(c(   seq(max(allParam[curIndex]-2*mult,  min.v),min(allParam[curIndex]+2*mult,  max.v),length.out=10)))



      if (method.rr=="COR")             {if (curParam == "R2" | curParam == "R1") mult2 <- 0.1
      if (curParam == "b")                     mult2  <- 2}
      if (method.rr=="COR") if(allParam[curIndex]  >= min.v & allParam[curIndex] <= max.v) poss.vals <- sort(c(allParam[curIndex],seq(max(allParam[curIndex]-mult2*mult,min.v),min(allParam[curIndex]+mul2*mult,max.v),length.out=10)))
      if (method.rr=="COR") if( (allParam[curIndex] < min.v  | allParam[curIndex] > max.v)) poss.vals <- sort(c(seq(max(allParam[curIndex]-mult2*mult,min.v),min(allParam[curIndex]+mul2*mult,max.v),length.out=10)))


      if (curIndex <   length(allParam) & method.rr=="TP") poss.vals <- sort(c(allParam[curIndex],seq(max(allParam[curIndex]-0.1*mult,0.000001),min(allParam[curIndex]+0.1*mult,0.9999),length.out=10)))
      if (curIndex ==  length(allParam) & method.rr=="TP") poss.vals <- sort(c(allParam[curIndex],seq    (allParam[curIndex]-0.1*mult,            allParam[curIndex]+0.1*mult,        length.out=10)))

      npv       <- length(poss.vals)
      func.vals <- matrix(nrow=npv,ncol=1)
      for (i in 1:npv) {allParam[curIndex]=poss.vals[i]
      if (method.rr=="OP" & altParam==0) func.vals[i] <- fL.OP(NA,allParam[1:s],allParam[(s+1):(2*s)],allParam[2*s+1],x1,x2,x3,n1,n2,n3)
      if (method.rr=="OP" & altParam==1) func.vals[i] <- fL.OP(allParam[1],NA,allParam[(2):(s+1)],allParam[s+2],x1,x2,x3,n1,n2,n3)
      if (method.rr=="CO")               func.vals[i] <- fL.CO(allParam[length(allParam)],allParam[1],allParam[-c(1,length(allParam))],Y,Z,tr)
      if (method.rr=="COR")              func.vals[i] <- fL.COR(allParam[length(allParam)],allParam[1],allParam[-c(1,length(allParam))],Y,Z,tr)
      if (method.rr=="TP")               func.vals[i] <- fL.TP(allParam[1:s],allParam[length(allParam)],allParam[(s+1):(2*s)],allParam[(2*s+1):(3*s)],x1,x2,x3,x4,n1,n2,n3,n4)
      }
      allParam[curIndex] <-poss.vals[which.min(func.vals)]
    }
    if (method.rr=="OP" & altParam==0) newMin <- fL.OP(NA,allParam[1:s],allParam[(s+1):(2*s)],allParam[2*s+1],x1,x2,x3,n1,n2,n3)
    if (method.rr=="OP" & altParam==1) newMin <- fL.OP(allParam[1],NA,  allParam[(2):(s+1)],  allParam[s+2],  x1,x2,x3,n1,n2,n3)
    if (method.rr=="CO")               newMin <- fL.CO(allParam[length(allParam)],allParam[1],allParam[-c(1,length(allParam))],Y,Z,tr)
    if (method.rr=="COR")              newMin <- fL.COR(allParam[length(allParam)],allParam[1],allParam[-c(1,length(allParam))],Y,Z,tr)
    if (method.rr=="TP")               newMin <- fL.TP(allParam[1:s],allParam[length(allParam)],allParam[(s+1):(2*s)],allParam[(2*s+1):(3*s)],x1,x2,x3,x4,n1,n2,n3,n4)
    if (newMin==prevMin) mult <- 0.9*mult
    prevMin <- newMin
    if (mult < 0.001) break
  }
  ok <- ifelse(mult < 0.1,1,0)
  if (mult < 0.1){
    if (method.rr=="OP" & altParam==0 & allRes==0 & fix.r==0) out <- allParam[2*s+1]
    if (method.rr=="OP" & altParam==0 & allRes==0 & fix.r==1) out <- list("pa"=allParam[1:s],"p3"=allParam[(s+1):(2*s)],"r"=allParam[2*s+1],"t"=NA)
    if (method.rr=="OP" & altParam==1 & allRes==0 & fix.r==0) out <- allParam[s+2]
    if (method.rr=="OP" & altParam==1 & allRes==0 & fix.r==1) out <- list("pa"=NA,           "p3"=allParam[(2):(s+1)],  "r"=allParam[s+2],  "t"=allParam[1])
    if (method.rr=="CO"               & allRes==0 & fix.r==0) out <- allParam[length(allParam)]
    if (method.rr=="CO"               & allRes==0 & fix.r==1) out <- list("r"=allParam[length(allParam)],"k"=allParam[1],"b"=allParam[-c(1,length(allParam))])
    if (method.rr=="COR"              & allRes==0 & fix.r==0) out <- allParam[length(allParam)]
    if (method.rr=="COR"              & allRes==0 & fix.r==1) out <- list("R1"=allParam[length(allParam)],"R2"=allParam[1],"b"=allParam[-c(1,length(allParam))])
    if (method.rr=="TP"               & allRes==0 & fix.r==0) out <- allParam[length(allParam)]
    if (method.rr=="TP"               & allRes==0 & fix.r==1) out <- list("tp"=allParam[1:s],"r"=allParam[length(allParam)],"p3"=allParam[(s+1):(2*s)],"p4"=allParam[(2*s+1):(3*s)])


    if (method.rr=="OP" & altParam==0 & allRes==1 & fix.r==1) out <- c(allParam[1:s],allParam[(s+1):(2*s)])
    if (method.rr=="OP" & altParam==0 & allRes==1 & fix.r==0) out <- c(allParam[1:s],allParam[(s+1):(2*s)],allParam[2*s+1])
    if (method.rr=="OP" & altParam==1 & allRes==1 & fix.r==1) out <- c(allParam[1],allParam[(2):(s+1)])
    if (method.rr=="OP" & altParam==1 & allRes==1 & fix.r==0) out <- c(allParam[1],allParam[(2):(s+1)],allParam[s+2])
    if (method.rr=="CO"               & allRes==1 & fix.r==1) out <- allParam[-length(allParam)]
    if (method.rr=="CO"               & allRes==1 & fix.r==0) out <- allParam
    if (method.rr=="COR"              & allRes==1 & fix.r==1) out <- allParam[-length(allParam)]
    if (method.rr=="COR"              & allRes==1 & fix.r==0) out <- allParam
    if (method.rr=="TP"               & allRes==1 & fix.r==1) out <- allParam[-length(allParam)]
    if (method.rr=="TP"               & allRes==1 & fix.r==0) out <- allParam

  }
  list("out"=out,"ok"=ok)
}



optim.ideal <- function(pa=NA,p3=NA,t=NA,r=NA,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="OP",x1,x2,x3,x4,n1,n2,n3,n4,n,Z,Y,tr,fix.r=0,allRes=0){
  s      <- length(x1)
  r.out  <- r
  t.out  <- NA
  pa.out <- rep(NA,s)
  p3.out <- rep(NA,s)
  k.out  <- NA
  b.out  <- rep(NA,length(b))
  p4.out <- rep(NA,s)
  tp.out <- rep(NA,s)
  EST    <- NA
  out    <- NA
  ll.end <- 99999999999
  ok <- 1

  if (method.rr=="OP" & altParam==0){

    if (fix.r==1){
      for (sss in 1:s){
        resb     <- nleqslv2(x=c(pa[sss],p3[sss]),fn=fL.pap3.parms,parms=c(1,x1[sss],x2[sss],x3[sss],n1[sss],n2[sss],n3[sss],r,1),positive=TRUE)
        ll.start <- fL.pap3.parms(c(pa[sss],p3[sss]),parms=c(1,x1[sss],x2[sss],x3[sss],n1[sss],n2[sss],n3[sss],r,0))
        try(ll.end   <- fL.pap3.parms(resb$root,         parms=c(1,x1[sss],x2[sss],x3[sss],n1[sss],n2[sss],n3[sss],r,0)),silent=TRUE)
        if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0) ok <- 0
        if (ok==1) pa.out[sss]  <- resb$root[1]
        if (ok==1) p3.out[sss]  <- resb$root[2]
        if (ok==1)             out <- list("pa"=pa.out,"p3"=p3.out)
        if (ok==1 & allRes==1) out <- c(pa.out,p3.out)
        if (ok==1) if (max(c(pa.out-r*p3.out/2,pa.out+r*p3.out/2,p3.out),na.rm=T) > 1 | min(c(pa.out-r*p3.out/2,pa.out+r*p3.out/2,p3.out),na.rm=T) < 0) ok <- 0
      }
    }

    if (fix.r==0){
      resb      <- nleqslv2(x=c(pa,p3,r+5),fn=fL.pap3r.parms.m5,parms=c(length(x1),x1,x2,x3,n1,n2,n3,1),positive=TRUE)
      ll.start  <- fL.pap3r.parms.m5(c(pa,p3,r+5),parms=c(length(x1),x1,x2,x3,n1,n2,n3,0))
      try(ll.end    <- fL.pap3r.parms.m5(resb$root,   parms=c(length(x1),x1,x2,x3,n1,n2,n3,0)),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0)ok <- 0
      if (ok==1)             r.out <- out  <- resb$root[2*s+1] -5
      if (allRes==1 & ok==1) out  <- resb$root - c(rep(0,2*s),5)
      if (ok == 1) {pa.out <- resb$root[1:s]
      p3.out <- resb$root[(s+1):(2*s)]}
      if (ok==1) if (max(c(pa.out-r.out*p3.out/2,pa.out+r.out*p3.out/2,p3.out),na.rm=T) > 1 | min(c(pa.out-r.out*p3.out/2,pa.out+r.out*p3.out/2,p3.out),na.rm=T) < 0) ok <- 0
    }

  }


  if (method.rr=="OP" & altParam==1){

    if (fix.r==1){
      resb     <- nleqslv2(x=c(t,p3),fn=fL.tp3.parms,parms=c(s,x1,x2,x3,n1,n2,n3,r,1),positive=TRUE)
      ll.start <- fL.tp3.parms(c(t,p3),  parms=c(s,x1,x2,x3,n1,n2,n3,r,0))
      try(ll.end   <- fL.tp3.parms(resb$root,parms=c(s,x1,x2,x3,n1,n2,n3,r,0)),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0) ok <- 0
      if (ok==1) t.out   <- resb$root[1]
      if (ok==1) p3.out  <- resb$root[-1]
      if (ok==1)             out <- list("t"=t.out,"p3"=p3.out)
      if (ok==1 & allRes==1) out <- c(t.out,p3.out)
      if (ok==1) if (max(c(t.out*p3.out - r*p3.out/2,t.out*p3.out + r*p3.out/2,p3.out),na.rm=T) > 1 | min(c(t.out*p3.out - r*p3.out/2,t.out*p3.out + r*p3.out/2,p3.out),na.rm=T) < 0) ok <- 0
    }

    if (fix.r==0){
      resb      <- nleqslv2(x=c(t,p3,r+5),fn=fL.tp3r.parms.m5,parms=c(length(x1),x1,x2,x3,n1,n2,n3,1),positive=TRUE)
      ll.start  <- fL.tp3r.parms.m5(c(t,p3,r+5),parms=c(length(x1),x1,x2,x3,n1,n2,n3,0))
      try(ll.end    <- fL.tp3r.parms.m5(resb$root,   parms=c(length(x1),x1,x2,x3,n1,n2,n3,0)),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0)ok <- 0
      if (ok==1)             out   <- resb$root[s+2]-5
      if (ok==1)             {t.out   <- resb$root[1]
      p3.out  <- resb$root[1+c(1:s)]
      r.out   <- out}
      if (allRes==1 & ok==1) out <- resb$root - c(rep(0,s+1),5)
      if (ok==1) if (max(c(t.out*p3.out - r.out*p3.out/2,t.out*p3.out + r.out*p3.out/2,p3.out),na.rm=T) > 1 | min(c(t.out*p3.out - r.out*p3.out/2,t.out*p3.out + r.out*p3.out/2,p3.out),na.rm=T) < 0) ok <- 0

    }

  }


  if (method.rr=="CO"){

    if (fix.r==1){
      resb     <- nleqslv2(x=c(k,b),fn=fL.kb.parms,parms=c(length(Y),r,1,Y,tr,c(Z)))
      ll.start <- fL.kb.parms(c(k,b),  parms=c(length(Y),r,0,Y,tr,c(Z)))
      try(ll.end   <- fL.kb.parms(resb$root,parms=c(length(Y),r,0,Y,tr,c(Z))),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0) ok <- 0
      if (ok==1) k.out   <- resb$root[1]
      if (ok==1) b.out  <- resb$root[-1]
      if (ok==1)             out <- list("k"=k.out,"b"=b.out)
      if (ok==1 & allRes==1) out <- c(k.out,b.out)
      if (ok==1) if (max(k.out-r/2,k.out+r/2,na.rm=T) > 1 | min(k.out-r/2,k.out+r/2,na.rm=T) < 0 ) ok <- 0
    }

    if (fix.r==0){
      resb      <- nleqslv2(x=c(k,b,r),fn=fL.kbr.parms,parms=c(length(Y),1,Y,tr,c(Z)))
      ll.start  <- fL.kbr.parms(c(k,b,r),parms=c(length(Y),0,Y,tr,c(Z)))
      try(ll.end    <- fL.kbr.parms(resb$root,   parms=c(length(Y),0,Y,tr,c(Z))),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0)ok <- 0
      if (ok==1)             out  <- resb$root[length(c(k,b))+1]
      if (ok==1)             {k.out  <- resb$root[1]
      r.out  <- out}
      if (allRes==1 & ok==1) out <- resb$root
      if (ok==1) if (max(k.out-r.out/2,k.out+r.out/2,na.rm=T) > 1 | min(k.out-r.out/2,k.out+r.out/2,na.rm=T) < 0 ) ok <- 0
    }

  }



  if (method.rr=="COR"){

    if (fix.r==1){
      resb     <- nleqslv2(x=c(R2,b),fn=fL.kb.parms2,parms=c(length(Y),R1,1,Y,tr,c(Z)))
      ll.start <- fL.kb.parms2(c(R2,b),  parms=c(length(Y),R1,0,Y,tr,c(Z)))
      try(ll.end   <- fL.kb.parms2(resb$root,parms=c(length(Y),R1,0,Y,tr,c(Z))),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0) ok <- 0
      if (ok==1) k.out   <- resb$root[1]
      if (ok==1) b.out  <- resb$root[-1]
      if (ok==1)             out <- list("R2"=k.out,"b"=b.out)
      if (ok==1 & allRes==1) out <- c(k.out,b.out)
    }

    if (fix.r==0){
      resb      <- nleqslv2(x=c(R2,b,R1),fn=fL.kbr.parms2,parms=c(length(Y),1,Y,tr,c(Z)))
      ll.start  <- fL.kbr.parms2(c(R2,b,R1),parms=c(length(Y),0,Y,tr,c(Z)))
      try(ll.end    <- fL.kbr.parms2(resb$root,   parms=c(length(Y),0,Y,tr,c(Z))),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0)ok <- 0
      if (ok==1)             out  <- resb$root[length(c(R2,b))+1]
      if (allRes==1 & ok==1) out <- resb$root
    }

  }


  if (method.rr=="TP"){

    if (fix.r==1){
      for (sss in 1:s){
        resb     <- nleqslv2(x=c(tp[sss],p3[sss],p4[sss]),fn=fL.tp01p02.parms,parms=c(1,x1[sss],x2[sss],x3[sss],x4[sss],n1[sss],n2[sss],n3[sss],n4[sss],r,1),positive=TRUE)
        ll.start <- fL.tp01p02.parms(c(tp[sss],p3[sss],p4[sss]),parms=c(1,x1[sss],x2[sss],x3[sss],x4[sss],n1[sss],n2[sss],n3[sss],n4[sss],r,0))
        try(ll.end   <- fL.tp01p02.parms(resb$root,         parms=c(1,x1[sss],x2[sss],x3[sss],x4[sss],n1[sss],n2[sss],n3[sss],n4[sss],r,0)),silent=TRUE)
        if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0) ok <- 0
        if (ok==1) tp.out[sss]  <- resb$root[1]
        if (ok==1) p3.out[sss]  <- resb$root[2]
        if (ok==1) p4.out[sss]  <- resb$root[3]
        if (ok==1)             out <- list("tp"=tp.out,"p3"=p3.out,"p4.out"=p4.out)
        if (ok==1 & allRes==1) out <- c(tp.out,p3.out,p4.out)
      }
    }

    if (fix.r==0){
      resb      <- nleqslv2(x=c(tp,p3,p4,r+5),fn=fL.tp01p02r.parms.m5,parms=c(length(x1),x1,x2,x3,x4,n1,n2,n3,n4,1),positive=TRUE)
      ll.start  <- fL.tp01p02r.parms.m5(c(tp,p3,p4,r+5),parms=c(length(x1),x1,x2,x3,x4,n1,n2,n3,n4,0))
      try(ll.end    <- fL.tp01p02r.parms.m5(resb$root,   parms=c(length(x1),x1,x2,x3,x4,n1,n2,n3,n4,0)),silent=TRUE)
      if ((sum(resb$f.root>0.01)>1) | sum(is.na(resb$f.root))>0 | (ll.end>ll.start) | sum(is.na(resb$root))>0)ok <- 0
      if (ok==1)             out  <- resb$root[3*s+1] -5
      if (allRes==1 & ok==1) out  <- resb$root - c(rep(0,3*s),5)
    }

  }


  list("out"=out,"ok"=ok)
}

optim.DSS <- function(pa=NA,p3=NA,t=NA,r=NA,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="OP",x1,x2,x3,x4,n1,n2,n3,n4,n,Z,Y,tr,fix.r=0,allRes=0,use.force=0){
  ideal = 1
  optim.out <- optim.ideal(pa=pa,p3=p3,t=t,r=r,k=k,b=b,p4=p4,tp=tp,R1=R1,R2=R2,altParam=altParam,method.rr=method.rr,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,n=n,Z=Z,Y=Y,tr=tr,fix.r=fix.r,allRes=allRes)
  if (optim.out$ok==0 | use.force==1) {ideal=0
  optim.out <- optim.force(pa=pa,p3=p3,t=t,r=r,k=k,b=b,p4=p4,tp=tp,R1=R1,R2=R2,altParam=altParam,method.rr=method.rr,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,n=n,Z=Z,Y=Y,tr=tr,fix.r=fix.r,allRes=allRes)}
  list("out"=optim.out$out,"ideal"=ideal)
}



####################################################################################################################
#####
##### Wrapper functions for finding the conditional MLE
#####
####################################################################################################################



opt.OP1 <- function(r,pa.start,p3.start,x1,x2,x3,n1,n2,n3,forceRoot=0,allRes=0,use.force=0){
  optim.DSS(pa=pa.start,p3=p3.start,t=NA,r=r,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="OP",x1=x1,x2=x2,x3=x3,x4=rep(0,length(x3)),n1=n1,n2=n2,n3=n3,n4=rep(0,length(x3)),n=NA,Z=NA,Y=NA,tr=NA,fix.r=1,allRes=allRes,use.force=use.force)$out
}


opt.OP2 <- function(r,t.start,p3.start,x1,x2,x3,n1,n2,n3,forceRoot=0,allRes=0,use.force=0){
  optim.DSS(pa=NA,p3=p3.start,t=t.start,r=r,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=1,method.rr="OP",x1=x1,x2=x2,x3=x3,x4=rep(0,length(x3)),n1=n1,n2=n2,n3=n3,n4=rep(0,length(x3)),n=NA,Z=NA,Y=NA,tr=NA,fix.r=1,allRes=allRes,use.force=use.force)$out
}


opt.TP <- function(r,t.start,p01.start,p02.start,x1,x2,x01,x02,n1,n2,n01,n02,allRes=0,use.force=0){
  optim.DSS(pa=NA,p3=p01.start,t=NA,r=r,k=NA,b=NA,p4=p02.start,tp=t.start,R1=NA,R2=NA,altParam=0,method.rr="TP",x1=x1,x2=x2,x3=x01,x4=x02,n1=n1,n2=n2,n3=n01,n4=n02,n=NA,Z=NA,Y=NA,tr=NA,fix.r=1,allRes=allRes,use.force=use.force)$out
}



opt.CO <- function(r,k.start,b.start,tr,Z,Y,allRes=0,use.force=0){
  optim.DSS(pa=NA,p3=NA,t=NA,r=r,k=k.start,b=b.start,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="CO",NA,NA,NA,NA,NA,NA,NA,NA,n=length(Y),Z=Z,Y=Y,tr=tr,fix.r=1,allRes=allRes,use.force=use.force)$out
}


opt.COR <- function(R1,R2.start,b.start,tr,Z,Y,allRes=0,use.force=0){
  optim.DSS(pa=NA,p3=NA,t=NA,r=NA,k=NA,b=b.start,p4=NA,tp=NA,R1=R1,R2=R2.start,altParam=0,method.rr="COR",NA,NA,NA,NA,NA,NA,NA,NA,n=length(Y),Z=Z,Y=Y,tr=tr,fix.r=1,allRes=allRes,use.force=use.force)$out
}


####################################################################################################################
#####
#####Wrapper function for getting the uncondition MLE of r
#####
####################################################################################################################


getEst.OP <- function(init.pa=NA,init.t=NA,init.p3,init.r,x1,x2,x3,n1,n2,n3,altParam=0,unCond=0,allRes=0,forceRoot=0,use.force=0){
  optim.DSS(pa=init.pa,p3=init.p3,t=init.t,r=init.r,k=NA,b=NA,p4=NA,tp=NA,R1=NA,R2=NA,altParam=altParam,method.rr="OP",x1=x1,x2=x2,x3=x3,x4=rep(0,length(x3)),n1=n1,n2=n2,n3=n3,n4=rep(0,length(x3)),n=NA,Z=NA,Y=NA,tr=NA,fix.r=0,allRes=allRes,use.force=use.force)$out
}

getEst.TP <- function(init.p4,init.tp,init.p3,init.r,x1,x2,x3,x4,n1,n2,n3,n4,allRes=0,use.force=0){
  optim.DSS(pa=NA,p3=init.p3,t=NA,r=init.r,k=NA,b=NA,p4=init.p4,tp=init.tp,R1=NA,R2=NA,altParam=0,method.rr="TP",x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,n=NA,Z=NA,Y=NA,tr=NA,fix.r=0,allRes=allRes,use.force=use.force)$out
}


getEst.CO <- function(init.k,init.b,init.r,Y,Z,tr,allRes=0,use.force=0){
  optim.DSS(pa=NA,p3=NA,t=NA,r=init.r,k=init.k,b=init.b,p4=NA,tp=NA,R1=NA,R2=NA,altParam=0,method.rr="CO",x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,n=length(Y),Z=Z,Y=Y,tr=tr,fix.r=0,allRes=allRes,use.force=use.force)$out
}

getEst.COR <- function(init.R2,init.b,init.R1,Y,Z,tr,allRes=0,use.force=0){
  optim.DSS(pa=NA,p3=NA,t=NA,r=NA,k=NA,b=init.b,p4=NA,tp=NA,R1=init.R1,R2=init.R2,altParam=0,method.rr="COR",x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,n=length(Y),Z=Z,Y=Y,tr=tr,fix.r=0,allRes=allRes,use.force=use.force)$out
}


####################################################################################################################
#####
#####Calculates the Z statistic for a given set of values
#####
####################################################################################################################



getZ2.OP1 <- function(pa,p3,r,x1,x2,x3,n1,n2,n3,unCond=0){
  I   <- fI.OP1(pa,p3,r,n1,n2,n3)
  dL  <- dL.OP1(pa,p3,r,x1,x2,x3,n1,n2,n3)
  Z   <- -9999999
  if (!is.na(I[1,1])){
    si  <- solve(I[-1,-1])
    if (unCond==0) v   <- I[1,1]-I[1,-1]%*%si%*%I[-1,1]
    if (unCond==1) v   <- I[1,1]
    if (v>0 & unCond==0)  Z   <- (dL[1]- I[1,-1]%*%si%*%dL[-1])/sqrt(v)
    if (v>0 & unCond==1)  Z   <- (dL[1])/sqrt(v)
  }
  Z
}

getZ2.OP2 <- function(t,p3,r,x1,x2,x3,n1,n2,n3,unCond=0){
  I   <- fI.OP2(t,p3,r,n1,n2,n3)
  dL  <- dL.OP2(t,p3,r,x1,x2,x3,n1,n2,n3)
  Z   <- -9999999
  if (!is.na(I[1,1])){
    si  <- solve(I[-1,-1])
    if (unCond==0) v   <- I[1,1]-I[1,-1]%*%si%*%I[-1,1]
    if (unCond==1) v   <- I[1,1]
    if (v>0 & unCond==0)  Z   <- (dL[1]- I[1,-1]%*%si%*%dL[-1])/sqrt(v)
    if (v>0 & unCond==1)  Z   <- (dL[1])/sqrt(v)
  }
  Z
}

getZ2.TP <- function(r,t,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02,unCond=0){
  I   <- fI.TP(t,r,p01,p02,n1,n2,n01,n02)
  dL  <- dL.TP(r,t,p01,p02,x1,x2,x01,x02,n1,n2,n01,n02)
  Z   <- -9999999
  if (!is.na(I[1,1])){
    si  <- solve(I[-1,-1])
    if (unCond==0) v   <- I[1,1]-I[1,-1]%*%si%*%I[-1,1]
    if (unCond==1) v   <- I[1,1]
    if (v>0 & unCond==0)  Z   <- (dL[1]- I[1,-1]%*%si%*%dL[-1])/sqrt(v)
    if (v>0 & unCond==1)  Z   <- (dL[1])/sqrt(v)
  }
  Z
}

getZ2.CO <- function(r,k,b,tr,Z,Y,unCond=0) {
  n <- length(Y)
  I   <- fI.CO(r,k,b,tr,n,Z)
  dL  <- dL.CO(r,k,b,tr,n,Z,Y)
  Z2   <- -9999999
  if (!is.na(I[1,1])){
    si  <- solve(I[-1,-1])
    if (unCond==0) v   <- I[1,1]-I[1,-1]%*%si%*%I[-1,1]
    if (unCond==1) v   <- I[1,1]
    if (v>0 & unCond==0)  Z2   <- (dL[1]- I[1,-1]%*%si%*%dL[-1])/sqrt(v)
    if (v>0 & unCond==1)  Z2   <- (dL[1])/sqrt(v)
  }
  Z2
}



getZ2.COR <- function(R1,R2,b,tr,Z,Y,unCond=0) {
  n <- length(Y)
  I   <- fI.COR(R1,R2,b,tr,n,Z)
  dL  <- dL.COR(R1,R2,b,tr,n,Z,Y)
  Z2   <- -9999999
  if (!is.na(I[1,1])){
    si  <- solve(I[-1,-1])
    if (unCond==0) v   <- I[1,1]-I[1,-1]%*%si%*%I[-1,1]
    if (unCond==1) v   <- I[1,1]
    if (v>0 & unCond==0)  Z2   <- (dL[1]- I[1,-1]%*%si%*%dL[-1])/sqrt(v)
    if (v>0 & unCond==1)  Z2   <- (dL[1])/sqrt(v)
  }
  Z2
}


####################################################################################################################
#####
#####Calculates Z values for a given value r (using MLE of other parameters)
#####
####################################################################################################################


getZ.OP1 <- function(r,x1,x2,x3,n1,n2,n3,unCond=0,forceRoot=0,sv=NA,use.force=0){
  if (sum(!is.na(sv))==0) sv  <- startVals.OP(x1,x2,x3,n1,n2,n3,r)
  sv2 <- opt.OP1(r,sv$pa,sv$p3,x1,x2,x3,n1,n2,n3,forceRoot=forceRoot,use.force=use.force)
  Z <- getZ2.OP1(sv2$pa,sv2$p3,r,x1,x2,x3,n1,n2,n3,unCond)
  list("Z"=Z,"sv2"=sv2)
}


getZ.OP2 <- function(r,x1,x2,x3,n1,n2,n3,unCond=0,forceRoot=0,sv=NA,use.force=0){
  if (sum(!is.na(sv))==0) sv  <- startVals.OP(x1,x2,x3,n1,n2,n3,r)
  sv2 <- opt.OP2(r,sv$t,sv$p3,x1,x2,x3,n1,n2,n3,forceRoot=forceRoot,use.force=use.force)
  Z <- getZ2.OP2(sv2$t,sv2$p3,r,x1,x2,x3,n1,n2,n3,unCond)
  list("Z"=Z,"sv2"=sv2)
}

getZ.TP <- function(r,x1,x2,x01,x02,n1,n2,n01,n02,unCond=0,sv=NA,use.force=0){
  if (sum(!is.na(sv))==0)  sv  <- startVals.TP(x1,x2,x01,x02,n1,n2,n01,n02,r)
  sv2 <- opt.TP(r,sv$t,sv$p01,sv$p02,x1,x2,x01,x02,n1,n2,n01,n02,use.force=use.force)
  Z <- getZ2.TP(r,sv2$t,sv2$p01,sv2$p02,x1,x2,x01,x02,n1,n2,n01,n02,unCond)
  list("Z"=Z,"sv2"=sv2)
}

getZ.CO <- function(r,tr,Z,Y,unCond=0,sv=NA,use.force=0){
  if (sum(!is.na(sv))==0) sv  <- startVals.CO(r,tr,Z,Y)
  sv2 <- opt.CO(r,sv$k,sv$b,tr,Z,Y,use.force=use.force)
  ZZ<- getZ2.CO(r,sv2$k,sv2$b,tr,Z,Y,unCond)
  list("Z"=ZZ,"sv2"=sv2)
}

getZ.COR <- function(R1,tr,Z,Y,unCond=0,sv=NA,use.force=0){
  if (sum(!is.na(sv))==0) sv  <- startVals.COR(R1,tr,Z,Y)
  sv2 <- opt.COR(R1,sv$R2,sv$b,tr,Z,Y,use.force=use.force)
  Z<- getZ2.COR(R1,sv2$R2,sv2$b,tr,Z,Y,unCond)
  list("Z"=Z,"sv2"=sv2)
}

getZ.OP <- function(r,x1,x2,x3,n1,n2,n3,altParam=0,unCond=0,forceRoot=0,sv=NA,use.force=0){
  if (altParam==0) {ZSV <- getZ.OP1(r,x1,x2,x3,n1,n2,n3,unCond,forceRoot=forceRoot,sv=sv,use.force=use.force)
  Z   <- ZSV$Z
  sv2 <- ZSV$sv2
  }
  if (altParam==1) {ZSV <- getZ.OP2(r,x1,x2,x3,n1,n2,n3,unCond,forceRoot=forceRoot,sv=sv,use.force=use.force)
  Z   <- ZSV$Z
  sv2 <- ZSV$sv2
  }
  list("Z"=Z,"sv2"=sv2)
}



####################################################################################################################
#####
#####Gets good starting values for the nuicance parametyers given a specified value of r
#####
####################################################################################################################


####Given a specific value of r, generate a set of starting values for the MLE for pa and p3
####
startVals2.OP <- function(x1,x2,x3,n1,n2,n3,r){
  if (x1+x2+x3==0) x3 <- x3+1
  d0 <- r*(x3/n3)
  L1 <- (x1*d0-n1-n2-2*x1)*d0+x1+x2
  L2 <- (n1+n2+x1)*d0-n1-n2-x1-x2
  L3 <- n1+n2
  L0 <- x1*d0*(1-d0)
  C  <- L2^3/(27*L3^3)-L1*L2/(6*L3^2)+L0/(2*L3)
  B  <- sign(C)*sqrt(L2^2/(9*L3^2)-L1/(3*L3))
  A  <- (1/3)*(pi+acos(C/B^3))
  p1.t <- 2*B*cos(A)-L2/(3*L3)
  p2.t <- p1.t+d0
  pa=(p1.t+p2.t)/2
  p3=x3/n3

  p3 <- ifelse(p3<0.000001,0.000001,p3)
  pa <- ifelse(pa < 0 | pa < r*p3/2 | is.na(pa), max(0,r*p3/2), pa)
  t  <- mean(pa/p3)
  t  <- ifelse(t > min((1-r*p3/2)/p3),min(1-r*p3/2)/p3,t)
  t  <- ifelse(t < r/2,r/2,t)
  list("pa"=pa,"p3"=p3,"t"=t)
}

startVals.OP <- function(x1,x2,x3,n1,n2,n3,r){
  s  <- length(x1)
  sv <- startVals2.OP(x1[1],x2[1],x3[1],n1[1],n2[1],n3[1],r)
  pa <- sv$pa
  p3 <- sv$p3
  if (s > 1) for (i in 2:s) { sv <- startVals2.OP(x1[i],x2[i],x3[i],n1[i],n2[i],n3[i],r)
  pa <- c(pa,sv$pa)
  p3 <- c(p3,sv$p3)}
  p3 <- ifelse(p3<0.000001,0.000001,p3)
  p1 <- pa - (r/2)*p3
  p2 <- pa + (r/2)*p3

  thresh.p <- 0.000001
  thresh.p2<- 0.99999

  p1 <- ifelse(p1 < thresh.p,thresh.p,p1)
  p2 <- ifelse(p1==thresh.p,p1 + r*p3,p2)
  p2 <- ifelse(p2 < thresh.p,thresh.p,p2)
  p1 <- ifelse(p2==thresh.p,p2 - r*p3,p1)

  p1 <- ifelse(p1 > thresh.p2,thresh.p2,p1)
  p2 <- ifelse(p1==thresh.p2, p1 + r*p3,p2)
  p2 <- ifelse(p2 > thresh.p2,thresh.p2,p2)
  p1 <- ifelse(p2==thresh.p2, p2 - r*p3,p1)
  pa <- (p1 + p2)/2
  pa <- ifelse(pa < 0 | pa < r*p3/2, max(0,r*p3/2), pa)


  N <- n1+n2+n3
  t <- sum(N*((pa-p3*r/2 + pa+p3*r/2 )/2))/sum(N*p3)
  t  <- ifelse(t > min((1-r*p3/2)/p3),min(1-r*p3/2)/p3,t)
  t  <- ifelse(t < r/2,r/2,t)

  list("pa"=pa,"p3"=p3,"t"=t)
}


startVals.TP <- function(x1,x2,x01,x02,n1,n2,n01,n02,r){
  sv <- startVals.OP(x1,x2,x3=(x01+x02),n1,n2,n3=(n01+n02),r)
  pcomb <- (x01+x02)/(n01+n02)
  p2 <- (sv$pa+r*sv$p3/2) * ((x02/n02)/pcomb)
  p1 <- (sv$pa-r*sv$p3/2) * ((x01/n01)/pcomb)
  p01 <- (x01/n01)
  p02 <- (x02/n02)
  t   <- p1/p01+p2/p02
  list("p01"=p01,"p02"=p02,"t"=t)
}

startVals.CO <- function(r,tr,Z,Y){
  sv <- -10
  if (sum(Y[tr==0])>0) sv <- log(sum(Y[tr==0])/sum(tr==0))
  b <- c(sv,rep(0,ncol(Z)-1))
  k <- (sum(Y[tr!=0])/sum(tr!=0))/(sum(Y[tr==0])/sum(tr==0))
  k <- ifelse(k < abs(r/2)+0.00001, abs(r/2)+0.00001, k)
  list("k"=k,"b"=b)
}


startVals.COR <- function(R1,tr,Z,Y){
  sv <- -10
  if (sum(Y[tr==0])>0) sv <- log(sum(Y[tr==0])/sum(tr==0))
  b <- c(sv,rep(0,ncol(Z)-1))
  R2 <- (sum(Y[tr==2])/sum(tr==2))/(sum(Y[tr==0])/sum(tr==0))
  list("R2"=R2,"b"=b)
}


####################################################################################################################
#####
#####Identifies good starting values of r
#####
####################################################################################################################


startR.OP <- function(x1,x2,x3,n1,n2,n3){
  ###Option 1: Taylor Series Expansion
  p1       <- sum(x1)/sum(n1) + 0.000001
  p2       <- sum(x2)/sum(n2) + 0.000001
  p3       <- sum(x3)/sum(n3) + 0.000001
  #v        <- (p1*(1-p1)/n1 + p2*(1-p2)/n2)/p3^2 + ((p1-p2)^2/(p3^4))*(p3*(1-p3)/n3)
  #v2       <- 0.25*(1/p3)^4*(p3*(1-p3)/n3)*p1*(1-p1)/n1 + 0.25*(1/p3)^4*(p3*(1-p3)/n3)*p2*(1-p2)/n2 + (1/p3)^6*((p1-p2)^2)*2*(p3*(1-p3)/n3)^2
  #v3       <- v + v2
  #rs       <- ((x2+0.5/n2)/(n2+0.5/n2) - (x1+0.5/n1)/(n1+0.5/n1))/((x3+0.5/n3)/(n3+0.5/n3))
  #out      <- sum(rs[x3>0]/v3[x3>0])/sum(1/v3[x3>0])
  out <- (p2-p1)/p3
  if (is.na(out)) out <- 0
  out
}

startR.TP <- function(x1,x2,x3,x4,n1,n2,n3,n4){
  ###Option 1: Taylor Series Expansion
  p1       <- sum(x1)/sum(n1) + 0.000001
  p2       <- sum(x2)/sum(n2) + 0.000001
  p3       <- sum(x3)/sum(n3) + 0.000001
  p4       <- sum(x4)/sum(n4) + 0.000001
  v        <- (p1*(1-p1)/n1)/p3^2 + (p2*(1-p2)/n2)/p4^2 + (p1^2/p3^4)*(p3*(1-p3)/n3) + (p2^2/p4^4)*(p4*(1-p4)/n4)
  rs       <- ((x2+0.5/n2)/(n2+0.5/n2))/((x4+0.5/n4)/(n4+0.5/n4)) - ((x1+0.5/n1)/(n1+0.5/n1))/((x3+0.5/n3)/(n3+0.5/n3))
  out      <- sum(rs[x3>0]/v[x3>0])/sum(1/v[x3>0])
  if (is.na(out)) out <- 0
  out
}

startR.CO <- function(tr,Z,Y){
  p1 <- mean(Y[tr==1],na.rm=T)
  p2 <- mean(Y[tr==2],na.rm=T)
  p3 <- mean(Y[tr==0],na.rm=T)
  out <- (p2-p1)/p3
  if (is.na(out) | p3==0) out <- 0
  out
}

startR.COR <- function(tr,Z,Y){
  p1 <- mean(Y[tr==1],na.rm=T)
  p3 <- mean(Y[tr==0],na.rm=T)
  out <- p1/p3
  if (is.na(out) | p3==0) out <- 0
  out
}


####################################################################################################################
#####
##### Identifies the bounds for the score-based confidence intervals
#####
####################################################################################################################


####Find the value of r that gives a z statistic that equals qnorm(alpha2) or qnorm(1-alpha2)
####
findR.int <- function(r,alpha2,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,altParam=0,method.rr="OP",unCond=0,sv=NA,use.force=0)    {if (method.rr=="OP") out1 <- getZ.OP(r,x1,x2,x3,n1,n2,n3,altParam,unCond,forceRoot=1,sv=sv,use.force=use.force)
if (method.rr=="TP") out1 <- getZ.TP(r,x1,x2,x3,x4,n1,n2,n3,n4,unCond,forceRoot=1,sv=sv,use.force=use.force)
if (method.rr=="CO") out1 <- getZ.CO(r,tr,Z,Y,unCond,sv=sv,use.force=use.force)
if (method.rr=="COR")out1 <- getZ.COR(r,tr,Z,Y,unCond,sv=sv,use.force=use.force)
out <- abs(    abs(out1$Z) - abs(qnorm(alpha2))  )
sv2 <- out1$sv2
list("out"=out,"sv2"=sv2)}



boundR2    <- function(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,R2=NA,alpha2=0.025,type="UB",altParam=0,method.rr="OP",unCond=0,fast=0)   {
  out  <- NA
  ntry <- 0
  use.force <- 0
  while (ntry <= 2 & is.na(out)) {
    ntry=ntry+1
    if (ntry==2) use.force <- 1

    if (type=="UB") min.r    <- EST+0.001
    if (type=="UB") max.r    <- EST+2
    if (type=="UB") start.r  <- min.r
    if (type=="UB") last.r   <- max.r
    if (type=="LB") min.r    <- EST-2
    if (type=="LB") max.r    <- EST-0.001
    if (type=="LB" & method.rr=="COR") min.r <- max(min.r,0.001)
    if (type=="LB" & method.rr=="COR") max.r <- max(max.r,0.001)
    if (type=="LB") start.r  <- max.r
    if (type=="LB") last.r   <- min.r
    mult <- 0.01
    if (abs(EST) <= 0.05 & fast==0) mult <- 0.001
    poss.r  <- seq(start.r,last.r,by=ifelse(type=="LB",-1,1)*mult)
    poss.d  <- matrix(NA,ncol=1,nrow=length(poss.r))
    sv      <- NA
    sv.list <- list()
    minReq2 <- ifelse(!is.na(x1[1])  & (sum(x1==0)>0 | sum(x2==0)>0) ,3, ifelse(is.na(x1[1])  & (sum(Y[tr==1])==0 | sum(Y[tr==2])==0),3,1))
    aDec <- 0
    for (i in 1:length(poss.r)) {
      aaa=try(out2 <- findR.int(r=poss.r[i],alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,sv=sv,use.force=use.force),silent=TRUE)
      if (class(aaa) != "try-error") poss.d[i] <- out2$out
      if (class(aaa) != "try-error") sv.list[[i]] <- sv        <- out2$sv2
      if (class(aaa) != "try-error") if (i > 1) if (poss.d[i] < poss.d[i-1] & !is.na(poss.d[i]) & !is.na(poss.d[i-1]) ) aDec <- 1
      if (i > 2) if ( ( (poss.d[i] > poss.d[i-1] & (poss.d[i-1] < minReq2 | poss.d[i-2] < minReq2) & aDec == 1 ) | is.na(poss.d[i]) | class(aaa) == "try-error") & (i > 3) ) sv <- sv.list[[i-3]]
      if (i > 2) if (   (poss.d[i] > poss.d[i-1] & (poss.d[i-1] < minReq2 | poss.d[i-2] < minReq2) & aDec == 1 ) | is.na(poss.d[i]) | class(aaa) == "try-error") break
    }
    reachMax <- ifelse(i==length(poss.r),1,0)
    poss.r   <- seq(poss.r[i-2],poss.r[i],by=ifelse(type=="LB",-1,1)*mult*0.1)
    if (i==3 & is.na(poss.d[2]))  poss.r   <- seq(EST+ifelse(type=="LB",-1,1)*0.0001,poss.r[i],by=ifelse(type=="LB",-1,1)*0.001)
    poss.d   <- matrix(ncol=1,nrow=length(poss.r))
    for (i in 1:length(poss.r)) {
      aaa=try(out2 <- findR.int(r=poss.r[i],alpha2=alpha2,x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,tr=tr,Z=Z,Y=Y,altParam=altParam,method.rr=method.rr,unCond=unCond,sv=sv,use.force=use.force),silent=TRUE)
      if (class(aaa) != "try-error") poss.d[i] <- out2$out
      if (class(aaa) != "try-error") sv        <- out2$sv2
      if (i > 1) if (poss.d[i] > poss.d[i-1] | is.na(poss.d[i]) | class(aaa) == "try-error") break
    }
    minReq <- ifelse(!is.na(x1[1]) & (sum(x1==0)>0 | sum(x2==0)>0) ,1, ifelse(is.na(x1[1])  & (sum(Y[tr==1])==0 | sum(Y[tr==2])==0),3,0.25))
    out <- ifelse(reachMax==1 | poss.d[i-1] > minReq,NA,poss.r[i-1])
  }
  out}


boundR    <- function(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,alpha2=0.025, type="UB",altParam=0,       method.rr="TP",     unCond=0, fast=0)  {
  out <- boundR2(EST,x1,x2,x3,x4,n1,n2,n3,n4,tr,Z,Y,alpha2=alpha2,type=type,altParam=altParam,method.rr=method.rr,unCond=unCond, fast=0)
  out
}

####################################################################################################################
#####
##### Identifies the bounds for the likelihood-based confidence intervals
#####
####################################################################################################################


boundR.LRT <- function(init.r=NA,x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=NA,Z=NA,tr=NA,init.pa=NA,init.p3=NA,init.p4=NA,init.t=NA,init.tp=NA,init.k=NA,init.b=NA,init.R1=NA,init.R2=NA,method.rr="OP",altParam=0,alpha2=alpha/2,type="LB"){


  comp        <- function(r=NA,x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=NA,Z=NA,tr=NA,s=NA,EST.1=NA,init.pa=NA,init.p3=NA,init.p4=NA,init.t=NA,init.tp=NA,init.k=NA,init.b=NA,R1=NA,init.R2=NA,method.rr="OP",altParam=0) {
    if (method.rr=="OP" & altParam==0) EST.2     <- opt.OP1(r,init.pa,init.p3,x1,x2,x3,n1,n2,n3,forceRoot=1,allRes=1,use.force=1)
    if (method.rr=="OP" & altParam==1) EST.2     <- opt.OP2(r,init.t,init.p3,x1,x2,x3,n1,n2,n3,forceRoot=1,allRes=1,use.force=1)
    if (method.rr=="TP")               EST.2     <- opt.TP(r,init.t,init.p3,init.p4,x1,x2,x3,x4,n1,n2,n3,n4,allRes=1,use.force=1)
    if (method.rr=="CO")               EST.2     <- opt.CO(r,init.k,init.b,tr,Z,Y,allRes=1,use.force=1)
    if (method.rr=="COR")              EST.2     <- opt.COR(R1,init.R2,init.b,tr,Z,Y,allRes=1,use.force=1)
    if (method.rr=="OP" & altParam==0) out       <- -2*(fL.pap3r(EST.1,x1,x2,x3,n1,n2,n3)-fL.pap3r(c(EST.2,r),x1,x2,x3,n1,n2,n3))
    if (method.rr=="OP" & altParam==1) out       <- -2*(fL.tp3r( EST.1,x1,x2,x3,n1,n2,n3)-fL.tp3r(c( EST.2,r),x1,x2,x3,n1,n2,n3))
    if (method.rr=="TP"              ) out       <- -2*(fL.tp01p02r(EST.1,x1,x2,x3,x4,n1,n2,n3,n4)-fL.tp01p02r(c(EST.2,r),x1,x2,x3,x4,n1,n2,n3,n4))
    if (method.rr=="CO"              ) out       <- -2*(fL.kbr(EST.1,Y,Z,tr)-fL.kbr(c(EST.2,r),Y,Z,tr))
    if (method.rr=="COR"             ) out       <- -2*(fL.kbr2(EST.1,Y,Z,tr)-fL.kbr2(c(EST.2,R1),Y,Z,tr))
    out
  }

  s           <- length(x1)
  if (method.rr=="OP")               EST.1     <- getEst.OP(init.pa,init.t,init.p3,init.r,x1,x2,x3,n1,n2,n3,altParam,unCond=0,allRes=1,use.force=1)
  if (method.rr=="TP")               EST.1     <- getEst.TP(init.p4,init.tp,init.p3,init.r,x1,x2,x3,x4,n1,n2,n3,n4,allRes=1,use.force=1)
  if (method.rr=="CO")               EST.1     <- getEst.CO(init.k,init.b,init.r,Y,Z,tr,allRes=1,use.force=1)
  if (method.rr=="COR")              EST.1     <- getEst.COR(init.R2,init.b,init.R1,Y,Z,tr,allRes=1,use.force=1)

  if (type=="UB" & method.rr != "COR") poss.r <- seq(init.r, init.r+2,by= 0.01)
  if (type=="LB" & method.rr != "COR") poss.r <- seq(init.r, init.r-2,by=-0.01)
  if (type=="UB" & method.rr == "COR") poss.r <- seq(init.R1, init.R1+2,by= 0.01)
  if (type=="LB" & method.rr == "COR") poss.r <- seq(init.R1, min(0.001,init.R1/10), by=-0.01)


  npr    <- length(poss.r)
  v      <- matrix(nrow=npr,ncol=1)
  for (i in 1:npr) {
    if (method.rr=="OP" & altParam==0) v[i] <- comp(r=poss.r[i],x1=x1,x2=x2,x3=x3,x4=NA,n1=n1,n2=n2,n3=n3,n4=NA,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=init.pa,init.p3=init.p3,init.p4=NA,init.t=NA,init.tp=NA,init.k=NA,init.b=NA,method.rr="OP",altParam=0)
    if (method.rr=="OP" & altParam==1) v[i] <- comp(r=poss.r[i],x1=x1,x2=x2,x3=x3,x4=NA,n1=n1,n2=n2,n3=n3,n4=NA,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=NA,init.p3=init.p3,init.p4=NA,init.t=init.t,init.tp=NA,init.k=NA,init.b=NA,method.rr="OP",altParam=1)
    if (method.rr=="TP"              ) v[i] <- comp(r=poss.r[i],x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=NA,init.p3=init.p3,init.p4=init.p4,init.t=NA,init.tp=init.tp,init.k=NA,init.b=NA,method.rr="TP",altParam=0)
    if (method.rr=="CO"              ) v[i] <- comp(r=poss.r[i],x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=Y,Z=Z,tr=tr,s=s,EST.1=EST.1,init.pa=NA,init.p3=NA,init.p4=NA,init.tp=NA,init.k=init.k,init.b=init.b,method.rr="CO",altParam=0)
    if (method.rr=="COR"             ) v[i] <- comp(r=NA,x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=Y,Z=Z,tr=tr,s=s,EST.1=EST.1,init.pa=NA,init.p3=NA,init.p4=NA,init.tp=NA,init.k=NA,init.b=init.b,init.R2=init.R2,R1=poss.r[i],method.rr="COR",altParam=0)
    if (abs(v[i]) > qchisq(1-2*alpha2,1)  ) break}
  if (i > 1) {if (type=="UB") poss.r2 <- seq(poss.r[i-1], poss.r[i],by= 0.001)
  if (type=="LB") poss.r2 <- seq(poss.r[i-1], poss.r[i],by= -0.001)
  npr2    <- length(poss.r2)
  v2      <- matrix(nrow=npr2,ncol=1)
  for (i in 1:npr2) {
    if (method.rr=="OP" & altParam==0) v2[i] <- comp(r=poss.r2[i],x1=x1,x2=x2,x3=x3,x4=NA,n1=n1,n2=n2,n3=n3,n4=NA,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=init.pa,init.p3=init.p3,init.p4=NA,init.t=NA,init.tp=NA,init.k=NA,init.b=NA,method.rr="OP",altParam=0)
    if (method.rr=="OP" & altParam==1) v2[i] <- comp(r=poss.r2[i],x1=x1,x2=x2,x3=x3,x4=NA,n1=n1,n2=n2,n3=n3,n4=NA,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=NA,init.p3=init.p3,init.p4=NA,init.t=init.t,init.tp=NA,init.k=NA,init.b=NA,method.rr="OP",altParam=1)
    if (method.rr=="TP"              ) v2[i] <- comp(r=poss.r2[i],x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=NA,init.p3=init.p3,init.p4=init.p4,init.t=NA,init.tp=init.tp,init.k=NA,init.b=NA,method.rr="TP",altParam=0)
    if (method.rr=="CO"              ) v2[i] <- comp(r=poss.r2[i],x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=Y,Z=Z,tr=tr,s=s,EST.1=EST.1,init.pa=NA,init.p3=NA,init.p4=NA,init.tp=NA,init.k=init.k,init.b=init.b,method.rr="CO",altParam=0)
    if (method.rr=="COR"             ) v2[i] <- comp(r=NA,x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=Y,Z=Z,tr=tr,s=s,EST.1=EST.1,init.pa=NA,init.p3=NA,init.p4=NA,init.tp=NA,init.k=NA,init.b=init.b,init.R2=init.R2,R1=poss.r2[i],method.rr="COR",altParam=0)
    if (abs(v2[i]) > qchisq(1-2*alpha2,1)  ) break}
  }


  # if (i > 1 & !is.na(EST) & abs(EST) < 0.01) {if (type=="UB") poss.r3 <- seq(poss.r2[i-1], poss.r2[i],by=  0.0001)
  #                                             if (type=="LB") poss.r3 <- seq(poss.r2[i-1], poss.r2[i],by= -0.0001)
  # npr3    <- length(poss.r3)
  # v3      <- matrix(nrow=npr3,ncol=1)
  # for (i in 1:npr3) {
  #   if (method.rr=="OP" & altParam==0) v2[i] <- comp(r=poss.r3[i],x1=x1,x2=x2,x3=x3,x4=NA,n1=n1,n2=n2,n3=n3,n4=NA,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=init.pa,init.p3=init.p3,init.p4=NA,init.t=NA,init.tp=NA,init.k=NA,init.b=NA,method.rr="OP",altParam=0)
  #   if (method.rr=="OP" & altParam==1) v2[i] <- comp(r=poss.r3[i],x1=x1,x2=x2,x3=x3,x4=NA,n1=n1,n2=n2,n3=n3,n4=NA,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=NA,init.p3=init.p3,init.p4=NA,init.t=init.t,init.tp=NA,init.k=NA,init.b=NA,method.rr="OP",altParam=1)
  #  if (method.rr=="TP"              ) v2[i] <- comp(r=poss.r3[i],x1=x1,x2=x2,x3=x3,x4=x4,n1=n1,n2=n2,n3=n3,n4=n4,Y=NA,Z=NA,tr=NA,s=s,EST.1=EST.1,init.pa=NA,init.p3=init.p3,init.p4=init.p4,init.t=NA,init.tp=init.tp,init.k=NA,init.b=NA,method.rr="TP",altParam=0)
  #   if (method.rr=="CO"              ) v2[i] <- comp(r=poss.r3[i],x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=Y,Z=Z,tr=tr,s=s,EST.1=EST.1,init.pa=NA,init.p3=NA,init.p4=NA,init.tp=NA,init.k=init.k,init.b=init.b,method.rr="CO",altParam=0)
  #  if (method.rr=="COR"             ) v2[i] <- comp(r=NA,x1=NA,x2=NA,x3=NA,x4=NA,n1=NA,n2=NA,n3=NA,n4=NA,Y=Y,Z=Z,tr=tr,s=s,EST.1=EST.1,init.pa=NA,init.p3=NA,init.p4=NA,init.tp=NA,init.k=NA,init.b=init.b,init.R2=init.R2,R1=poss.r3[i],method.rr="COR",altParam=0)
  #  if (abs(v2[i]) > qchisq(1-2*alpha2,1)  ) break}
  # poss.r2=poss.r3
  # }


  if ( (poss.r2[i]==init.r+2 | poss.r2[i]==init.r-2)   & method.rr != "COR") poss.r2[i] <- NA
  if ( (poss.r2[i]==init.R1+2) & method.rr == "COR") poss.r2[i] <- NA
  poss.r2[i]
}

#######################################################################################################################################################################
###
#######################################################################################################################################################################


