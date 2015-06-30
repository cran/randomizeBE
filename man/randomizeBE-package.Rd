\name{randomizeBE-package}
\alias{randomizeBE-package}
\alias{randomizeBE}
\docType{package}
\title{
Create a randomisation list for cross-over studies.
  }
\description{
Main function \code{RL4()} randomizes subjects, patients in groups of 
(treatment) sequences.\cr
If a blocksize is given, the randomization will be done within blocks.\cr
If more then one blocksizes are given, the size of blocks will be chosen randomly
from the blocksizes.\cr
The randomization may be controlled by a Wald-Wolfowitz runs test.
Functions to obtain the p-value of that test are included.\cr\cr
Contains further two helper functions \code{sequences()} and \code{williams()} 
to get the sequences of commonly used designs in BE studies.
  }
\details{
\tabular{ll}{
Package: \tab randomizeBE\cr
Type: \tab Package\cr
Version: \tab 0.3-1\cr
Date: \tab 2012-09-21\cr
License: \tab GPL (>=2)\cr
  }
}
\author{
D. Labes\cr
with part of the code for using the runs test for randomness \cr
according to a suggestion by Helmut Schuetz.

Maintainer: D.Labes <detlewlabes@gmx.de>
}

