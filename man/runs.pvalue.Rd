\name{runs.pvalue}
\alias{runs.pvalue}
\title{
2-sided p-value of the runs test
}
\description{
The function calculates the 2-sided p-value of the Wald-Wolfowitz runs test
after dichotomizing the input vector
}
\usage{
runs.pvalue(y, pmethod = c("exact", "normal", "cc"))
}
\arguments{
  \item{y}{
Numeric vector of data values.
}
  \item{pmethod}{
A character string describing the method for the p-value calculation of the runs
test.\cr\cr
\code{"exact"} chooses the calculation via exact distribution of the # of runs.\cr
\code{"normal"} denotes the normal approximation like the function(s) \code{runs.test()}
of the packages \code{tseries} or \code{lawstat}.\cr
\code{"cc"} chooses the continuity correction to the large sample approximation 
like in the statistical software SPSS.
}
}
\details{
This function calculates the 2-sided p-value of the runs.test.\cr
The large sample approximations are an adaption from the codes for \code{runs.test()} 
found in the R-packages lawstat and tseries.\cr 
The aim of this own was to avoid the heavy footprint of both packages for this small
package.\cr\cr
The user can choose the application of a continuity correction to the normal
approximation like a SAS implementation \url{http://support.sas.com/kb/33/092.html} uses  
or like SPSS if n<50.\cr\cr
The exact distribution of runs and the p-value based on it are described in the manual
of SPSS "Exact tests" to be found f.i. \url{http://www.sussex.ac.uk/its/pdfs/SPSS_Exact_Tests_21.pdf}.\cr
If \code{pmethod="exact"} is chosen and n>30 and n1>12 and n2>12 (see \code{\link{pruns.exact}})
the continuity corrected version of the normal approximation is used to save time and memory.
}
\value{
Numeric p-value of the 2-sided test.
}
\author{
D. Labes\cr
adapted from \code{runs.test()} package lawstat\cr
Authors: Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao\cr
and from \code{runs.test()} package tseries\cr
Author: A. Trapletti
}
\seealso{
\code{\link{pruns.exact}}
}
\examples{
# alternating sequence 1,2,1,2 ...
# maybe seen as numeric representation of 'TR','RT' ...
# and is used in that way here in this package
x <- rep(c(1, 2), 6)
runs.pvalue(x, pmethod="normal")
# should give 0.002464631
# exact p-value
runs.pvalue(x, pmethod="exact")
# should give 0.004329004
# 
# same for 3 numbers (numeric representation of 3 sequences)
x <- rep(c(1, 2, 3),4)
runs.pvalue(x, pmethod="normal")
# should give 0.2502128
# i.e. is seen as compatible with a random sequence!
# exact p-value, default i.e. must not given exolicitely
runs.pvalue(x)
# should give 0.3212121
# i.e. is seen even more as compatible with a random sequence!
}
