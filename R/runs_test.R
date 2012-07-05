###############################################################################
# p.value of the 2-sided runs test 
# 
# Author: dlabes
###############################################################################
# adapted from runs.test package lawstat
# Authors: Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
# and runs.test package tseries
# Author(s): A. Trapletti
# only the 2-sided p-value is of concern for our purposes
# in package randomizeBE

runs.pvalue <- function(y){
  y <- na.omit(y)
  med <- median(y, na.rm = TRUE)
  # for values == median the preceeding value is used
  for (k in 2:length(y)) {
    if ((y[k] == med) & (y[k - 1] < med)) {
      y[k] = y[k - 1]
    }
    else if ((y[k] == med) & (y[k - 1] > med)) {
      y[k] = y[k - 1]
    }
  }
  # -1 and +1
  s  <- sign(y - med)
  n  <- length(s)
  # runs: see runs.test from package tseries
  R  <- 1 + sum(as.numeric(s[-1] != s[-n]))
  # number of above/ below
  n1 <- sum(s == +1)
  n2 <- sum(s == -1)
  E  <- 1 + 2*n1*n2/(n1 + n2)
  s2 <- (2*n1*n2*(2*n1*n2 - n1 - n2))/((n1 + n2)^2 * (n1 + n2 - 1))
  statistic <- (R - E)/sqrt(s2)
  pvalue    <- 2 * pnorm(-abs(statistic))
  pvalue
}
