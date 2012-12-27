\name{sequences}
\alias{sequences}
\title{
Obtain sequences for BE study designs
}
\description{
The function is intended as helper function to get the sequences for commonly 
used designs in BE studies.\cr 
These sequences may then subsequently used in \code{RL4()}.
}
\usage{
sequences(design, tmts = NULL)
}
\arguments{
  \item{design}{
A character value describing the study design. For crossover studies it is
"tmts x sequences x periods" without space.\cr
F.i. "2x2x2" is the standard 2-treatments-2-sequence-2-period crossover.\cr
The designs \code{"2x2x2", "3x3x3", "4x4x4"} may be abbreviated to 
\code{"2x2", "3x3", "4x4"}.
}
  \item{tmts}{
A character vector with the codes of the treatments.
}
}
\details{
This function was implemented because I couldn't remember f.i. all the six 
sequences of a "3x6x3" design if had to be coded by T1,T2 and R.\cr
The sequences are primarily build within the ABC... notation.\cr
If treatment codes are given the "AB..." in the sequences are replaced by these codes.\cr\cr
The following designs are implemented:\cr
\code{"parallel"} = 2-group parallel design\cr
\code{"2x2"} = classical 2-treatments-2-sequence-2-period crossover\cr
\code{"3x3"} = 3-treatments-3-sequence-3-period crossover (Latin square)\cr
\code{"3x6x3"} = 3-treatments-6-sequence-3-period crossover (Williams design)\cr
\code{"4x4"} = 4-treatments-4-sequence-4-period crossover (Latin square)\cr
\code{"2x2x3"} = 2-sequence-3-period replicate crossover\cr
\code{"2x2x4"} = 2-sequence-4-period full replicate crossover\cr
\code{"2x4x4"} = 4-sequence-4-period full replicate crossover\cr
\code{"2x3x3"} = 3-sequence-3-period partial (reference) replicate crossover\cr
\code{"2x4x2"} = 4-sequence-2-period replicate crossover, Baalams design.\cr\cr
The sequences for the \code{"3x3"} and \code{"4x4"} designs are randomly derived 
from the standard Latin squares (in ABC notation)\cr
\code{ABC, BCA, CAB}\cr
or in case of the \code{"4x4"} design from \cr
\code{ABCD, BDAC, CADB, DCBA}.
}
\value{
Returns a character vector with the sequences.
}
\references{
Jones B, Kenward MG (2003). 
"Design and Analysis of Cross-Over Trials" 2nd edition.
Chapman & Hall, London.
}
\author{
Detlew Labes
}
\note{
For the higher order designs (designs with more than 2 treatments or replicate
crossover designs, respectively) only a selection of possible designs are covered.\cr
The design sequences for \code{"4x4"} are not guaranteed to be a Williams design.
If need a Williams design use function \code{williams()} instead. 
}
\seealso{
\code{\link{williams}}, \code{\link{RL4}}
}
\examples{
# classical 2x2 crossover in TR notation (simple enough to remember
sequences("2x2", tmts=c("T","R"))
# 3-treatment-6-sequence-3-period Williams design in ABC notation
sequences("3x6x3")
# 3-treatment-3-sequence-3-period design with one Test and two References
sequences("3x3", tmts=c("T","R1","R2"))
# 4-treatment-4-period in TxRy notation, two Test and two Reference
sequences("4x4", tmts=c("T1","T2","R1","R2"))
}
