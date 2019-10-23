% --- Source file: options.Rd ---
\name{options}
\alias{options}
\docType{data}
\title{
 Options
}
\description{
 The list of options
}
\format{
      The format is: List of 8
       \describe{
           \item{alpha}{The two-sided error rate for the confidence interval. The default is 0.05.}
           \item{altParam}{A binary variable (i.e. 0 or 1) indivating whether to use the 
                           full likelihood (\code{altParam=0})
                           or assume the log-linear model holds true (\code{altParam=1}). 
                           The default is \code{altParam = 0}.
                           We note, \code{altParam=1} is equivalent to using \code{\link{DRRCI_adj}}  
                           with covars set to binary indicators for strata.}
           \item{deltaMethod}{A binary variable (i.e. 0 or 1) indicating whether the delta method 
                              should (\code{deltaMethod=1}) or should not (\code{deltaMethod=0}) be used. 
                              The default is \code{deltaMethod = 0}.}
           \item{estRR}{A binary variable (i.e. 0 or 1) indicating whether we should calculate the 
                        confidence intervals for each relative risk (\code{estRR=1}) or the difference 
                        in relative risks (\code{estRR=0}). The default is \code{estRR = 0}.}
           \item{fast}{A binary variable (i.e. 0 or 1)  indicating that the optimizations should be performed 
                       over a coarser grid to speed up the analysis (\code{fast=1}). 
                       Use with caution. The default is \code{fast = 0}.}
           \item{fast.c}{A binary variable (i.e. 0 or 1)  indicating that the c code should be called
                         instead of the R code for faster computing times. 
                         The default is \code{fast.c = 1}.}
           \item{LRT}{A binary variable (i.e. 0 or 1) indicating the type of confidence interval desired:
                      0 = a score-based confidence interval and 1 = a likelihood-based confidence interval. 
                      The default is \code{LRT = 0}.}
           \item{refPop}{NULL or a string to denote the target population when standardization is used.
                        The choices are 
                         \tabular{ll}{
                           Choice  \tab Description \cr
                           "All"   \tab All subjects \cr
                           "Cont"  \tab The control population \cr
                           "Ind"   \tab Indirect standardization for the delta method \cr
                           "Opt"   \tab An optimally chosen population \cr
                           "Treat" \tab The treated population
                         } 
                 The default is NULL unless \code{deltaMethod=1}, in that case, the default is "All".
           }
     }
 }
       
\details{Not all \code{options} are for all functions.
 Below is a table showing which \code{options} can be used.
                        \tabular{lccc}{
                           Option      \tab DRRCI_adj \tab DRRCI_counts           \tab DRRCI_unadj    \cr
                           alpha       \tab Yes       \tab Yes                    \tab Yes            \cr
                           altParam    \tab No        \tab Yes (3 groups)         \tab Yes (3 groups) \cr
                           deltaMethod \tab No        \tab Yes (3 groups)         \tab Yes (3 groups) \cr 
                           estRR       \tab Yes       \tab No                     \tab No             \cr
                           fast        \tab Yes       \tab Yes                    \tab Yes            \cr
                           fast.c      \tab Yes       \tab Yes                    \tab Yes            \cr
                           LRT         \tab Yes       \tab Yes                    \tab Yes            \cr
                           refPop      \tab No        \tab Yes                    \tab Yes      
                         } 

\bold{Notes:} \cr
1. The option \code{fast} is only valid when \code{LRT = 0}. \cr
2. With four treatment groups, \code{refPop = "Opt"} is not valid and will default to \code{refPop = "All"}.

}
%\examples{}
\keyword{misc}
\seealso{ \code{\link{DRRCI_adj}}, \code{\link{DRRCI_unadj}}, \code{\link{DRRCI_counts}} }

