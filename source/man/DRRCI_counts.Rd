\name{DRRCI_counts}
\alias{DRRCI_counts}
\title{Difference in Relative Risks}
\description{
This function estimates and provides a confidence interval for the difference between
two relative risks with count data. 
}

\usage{
DRRCI_counts(x1, n1, x2, n2, x0, n0, x0b=NULL, n0b=NULL, options=NULL)
}
\arguments{
\item{x1}{The number of events in the first treatment group. This can also be a vector
          if there are strata.}
\item{n1}{The number of subjects in the first treatment group. This can also be a vector
          if there are strata.}
\item{x2}{The number of events in the second treatment group. This can also be a vector
          if there are strata.}
\item{n2}{The number of subjects in the second treatment group. This can also be a vector
          if there are strata.}
\item{x0}{The number of events in the (first) control group. This can also be a vector
          if there are strata.}
\item{n0}{The number of subjects in the (first) control group. This can also be a vector
          if there are strata.}
\item{x0b}{NULL or the number of events in the (second) control group. This can also be a vector
          if there are strata.}
\item{n0b}{NULL or the number of subjects in the (second) control group. This can also be a vector
          if there are strata.}

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
 This function estimates and provides a confidence interval for the difference between
 two relative risks. If \code{x0b} is \code{NULL}, then there are three populations,
 with one control population, and 
 the goal to estimate the difference: \code{p1/p0 - p2/p0}, where \code{p0}, \code{p1}, and
 \code{p2} are the risks in each population.
 When \code{x0b} is not \code{NULL}, then there are four populations (two control populations),
 with \code{x0, n0} as the control population for the first treatment group and \code{x0b, n0b}
 as the control population for the second treatment group.
 With two control populations,  the goal is to estimate
 the difference: \code{p2/p0b - p1/p0}.
 
 If any group has zero events, then for the calculation of the lower and upper bounds with
 \code{options$deltaMethod = 1}, the zero will be replaced with 0.5.

}
\examples{
  # No strata
  DRRCI_counts(15, 29, 16, 40, 17, 31)
  DRRCI_counts(15, 29, 16, 40, 17, 31, 
               options=list(alpha=0.05,altParam=0))

  # With 3 strata
  DRRCI_counts(c(6, 4, 5), c(10, 9, 10), 
               c(5, 7, 4), c(13, 15, 12), 
               c(5, 7, 5), c(10, 10, 11))
  DRRCI_counts(c(6, 4, 5), c(10, 9, 10), 
               c(5, 7, 4), c(13, 15, 12), 
               c(5, 7, 5), c(10, 10, 11),
               options=list(alpha=0.05,altParam=0))


}
\seealso{ \code{\link{DRRCI_adj}}, \code{\link{DRRCI_unadj}} }

