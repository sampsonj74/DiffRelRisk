\name{DRRCI_unadj}
\alias{DRRCI_unadj}
\title{Difference in Relative Risks}
\description{
This function estimates and provides a confidence interval for the difference between
two relative risks with no covariates. 
}

\usage{
DRRCI_unadj(trtGroup, outcome, strata=NULL, options=NULL)
}
\arguments{
\item{trtGroup}{A vector of (3 or 4) treatment groups coded 0-2, or 1-4 (see details).}
\item{outcome}{A binary vector of outcomes coded as 0 and 1.}
\item{strata}{NULL or a vector of strata.}
\item{options}{List of options (see \code{\link{options}}).}
}
\value{
A list containing: 
 \tabular{ll}{
  Name \tab Desription \cr
  \code{EST}  \tab The estimated difference in relative risks \cr
  \code{LB}   \tab The lower bound for \code{EST} \cr
  \code{UB}   \tab The upper bound for \code{EST} \cr
 }
}
\details{
This function is essentially \code{\link{DRRCI_counts}}, but with a different interface.
If \code{trt} is coded 0-2, then it is assumed that 0 denotes the control group,
 1 denotes treatment group 1, and 2 denotes treatment group 2. If there are four groups
 (\code{trt} is coded 1-4), then \code{trt=3} denotes the control group for
 treatment group 1 (\code{trt=1}), and \code{trt=4} denotes the control group for
 treatment group 2 (\code{trt=2}).
}
\examples{
  set.seed(1234)
  n   <- 100
  trt <- sample(0:2, n, replace=TRUE)
  Y   <- sample(0:1, n, replace=TRUE)
  S   <- sample(1:3, n, replace=TRUE)

  DRRCI_unadj(trt, Y, options=list(deltaMethod=0, LRT=0))
  DRRCI_unadj(trt, Y, options=list(deltaMethod=0, LRT=1))
  DRRCI_unadj(trt, Y, options=list(deltaMethod=1, LRT=0))

  DRRCI_unadj(trt, Y, strata=S)
  DRRCI_unadj(trt, Y, strata=S, options=list(deltaMethod=0, LRT=0, alpha=0.05))
}
\seealso{ \code{\link{DRRCI_adj}}, \code{\link{DRRCI_counts}} }


