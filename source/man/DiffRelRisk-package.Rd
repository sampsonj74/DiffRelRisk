
\name{DiffRelRisk-package}
\alias{DiffRelRisk}
\docType{package}
\title{
 Confidence Intervals for the Difference Between Two Relative Risks
}
\description{
  Estimates and provides a confidence interval for the difference between two relative risks
}
\details{
 This package contains methods to estimate the confidence interval for the 
difference between two relative risks. 
Letting \code{p0}, \code{p1}, and \code{p2} be the probabilities of an event 
in three groups (i.e. control, treatment 1, treatment 2), these methods estimate 
a confidence interval for \code{r = p1/p0 - p2/p0}. 
These methods can handle small sample sizes, stratified studies, 
and individual-level covariates. 
They were developed specifically for vaccine trials to estimate the
 difference between two vaccine efficacies \code{(VE)},
 where \code{VE1 = 1- p1/p0}, \code{VE2 = 1- p2/p0} and \code{r = VE2 - VE1}. 
}
\references{
 
}
\author{Joshua Neil Sampson <joshua.sampson@nih.gov> and
        Mitchell H. Gail <gailm@exchange.nih.gov>}
\keyword{package, confidence interval, vaccine efficacy, relative risks, 
difference in relative risks, small samples}