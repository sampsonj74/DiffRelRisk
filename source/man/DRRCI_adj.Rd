\name{DRRCI_adj}
\alias{DRRCI_adj}
\title{Difference in Relative Risks}
\description{
This function estimates and provides a confidence interval for the difference between
two relative risks when adjusting for covariates. 
}

\usage{
DRRCI_adj(trtGroup, outcome, covars, options=NULL)
}
\arguments{
\item{trtGroup}{A vector of 3 treatment groups coded 0-2, with 0=control, 
                1=treatment group 1, and 2=treatment group 2.}
\item{outcome}{A binary vector of outcomes coded as 0 and 1.}
\item{covars}{A matrix of covariates}
\item{options}{List of options (see \code{\link{options}}).}
}
\value{
A list containing: 
 \tabular{ll}{
  Name \tab Desription \cr
  \code{EST}  \tab The estimated difference in relative risks \cr
  \code{LB}   \tab The lower bound for \code{EST} \cr
  \code{UB}   \tab The upper bound for \code{EST} \cr
  \code{beta} \tab The MLE (and asymptotic SE) for the coefficients for individual level covariates
 }
}
\details{
The objects \code{trtGroup}, \code{outcome}, and \code{covars} must all have the same 
number of subjects. Subjects with missing values will be removed.
}

\examples{
  set.seed(1234)
  n   <- 100
  trt <- sample(0:2, n, replace=TRUE)
  Y   <- sample(0:1, n, replace=TRUE)
  X   <- runif(n)

  DRRCI_adj(trt, Y, covars=X)
  DRRCI_adj(trt, Y, covars=X, options=list(alpha=0.05,LRT=0))

}
\seealso{ \code{\link{DRRCI_unadj}}, \code{\link{DRRCI_counts}} }

