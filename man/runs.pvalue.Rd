\name{runs.pvalue}
\alias{runs.pvalue}
\title{
2-sided p-value of the runs test
}
\description{
The function calculates the 2-sided p-value of the Wald-Wolfowitz runs test
after dichotomising the input vector
}
\usage{
runs.pvalue(y)
}
\arguments{
  \item{y}{
numeric vector of data values.
}
}
\details{
This is an adaption from the codes for \code{runs.test()} found in the R-packages  
lawstat and tseries to calculate the 2-sided p-value of the runs.test.\cr
The aim was to avoid the heavy footprint of both packages for this small
package.
}
\value{
numeric p-value of the 2-sided test.
}
\author{
D. Labes\cr
adapted from \code{runs.test()} package lawstat\cr
Authors: Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao\cr
and from \code{runs.test()} package tseries\cr
Author: A. Trapletti
}
\seealso{
\code{\link[tseries]{runs.test} tseries}, 
\code{\link[lawstat]{runs.test} lawstat}
}
\examples{
# alternating sequence 1,2,1,2 ...
# maybe seen as numeric representation of 'TR','RT' ...
# and is used in that way here in this package
x <- rep(c(1, 2), 6)
runs.pvalue(x)
# should give 0.002464631

# same for 3 numbers (numeric representation of 3 sequences)
x <- rep(c(1, 2, 3),4)
runs.pvalue(x)
# should give 0.2502128
# i.e. is seen as compatible with a random sequence!
}
