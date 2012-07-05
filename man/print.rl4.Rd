\name{print.rl4}
\alias{print.rl4}
\title{
S3 method \code{print} for class "rl4"
}
\description{
Prints the randomisation list including a summary if needed.
}
\usage{
\method{print}{rl4}(x, sumry=FALSE, ...)
}
\arguments{
  \item{x}{
Object of class "rl4".
}
  \item{sumry}{
If set to TRUE a summary of the randomisation will be shown.\cr
This summary contains the p-value of runs test of randomness if 2 sequences
are used.
}
  \item{\dots}{
Additional arguments. None used.
}
}
\value{
Returns invisible x.
}
\author{
D. Labes
}
\seealso{
\code{\link{RL4}}
}
\examples{
# block randomisation of 12 subjects for a reference replicate study
# and demonstration of the S3 print method
rl <- RL4(nsubj=12, blocksize=6, seqs=c("TRR","RTR","RRT"), seed=0)
print(rl, sumry=TRUE)
}

