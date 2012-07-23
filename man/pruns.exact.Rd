\name{pruns.exact}
\alias{pruns.exact}
\title{
Exact cumulative distribution function of runs test
}
\description{
This function calculates the exact cumulative conditional distribution of the
Wald-Wolfowitz runs. 
}
\usage{
pruns.exact(r, n1, n2, tail = c("2-sided", "lower", "upper"))
}
\arguments{
  \item{r}{
Number of runs observed.
}
  \item{n1}{
Number of +1 items in the sequence.
}
  \item{n2}{
Number of -1 items in the sequence.
}
  \item{tail}{
Tail of the cumulative distribution function.\cr
Default is the 2-tailed value.
}
}
\value{
Numeric value of the cumulative distribution function according to the chosen tail.
}
\author{
D. Labes
}
\note{
The 2-sided exact p-value is defined as P(abs(R-E(R)>=abs(r-E(R)).\cr
The lower (left) tail p-value is defined as P(R<=r).\cr
The upper (right) tail p-value is defined as P(R>=r).\cr
r is the observed value of the random variable R.

}
\seealso{
\code{\link{help}}
}
\examples{
# SPSS "Exact Tests": small sample example, exact p:  0.071
# x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)  
pruns.exact(r=3, n1=4, n2=6)
# 0.07142857
# left tail P(R<=3)=0.04761905
pruns.exact(r=3, n1=4, n2=6, tail="lower")
# right tail P(R>=3)=0.9904762
pruns.exact(r=3, n1=4, n2=6, tail="upper")
# or via runs.pvalue (2-sided)
x <- c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1)  
runs.pvalue(x, pmethod="ex")
}

