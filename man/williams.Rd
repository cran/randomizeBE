\name{williams}
\alias{williams}
\title{
Construct sequences of a Williams design
}
\description{
The function constructs the sequences of a Williams design via the algorithm
given by Sheehe and Bross.
}
\usage{
williams(ntmt = 4, tmts = NULL)
}
\arguments{
  \item{ntmt}{
Number of treatments. Must be an integer >1.
}
  \item{tmts}{
NULL or a vector of treatment codes.
}
}
\details{
A Williams design is a crossover design in which each subject receives each 
treatment. The design is balanced over periods. The design has additionally the 
further feature that every treatment follows every other treatment the same 
number of times. It is said that these designs are balanced for first order 
carry over effects.\cr
For an even number of treatments the design is a Latin square. For an odd number
the design is a combination of two Latin squares.\cr\cr
Although the balance to first order carry over effects is not absolutely necessary
in well planned BE studies with sufficient washout the Williams designs were 
traditionally used in \code{"3x6x3"} \code{"4x4"} crossover studies.\cr\cr
The sequences are created originally within the ABC... notation. If \code{tmts}
are given the sequences are returned based on these treatment codes. 
}
\value{
Returns a character vector of the sequences.
}
\references{
Williams, E. J. (1949)
"Experimental designs balanced for the estimation of residual 
effects of treatments"
Australian J. of Scientific Research, Ser. A 2, 149-168.

Sheehe PR, Bross IDJ (1961)
"Latin Squares to Balance Immediate Residual and Other Effects." 
Biometrics, 17, 405-414.

Jones B, Kenward MG (2003). 
"Design and Analysis of Cross-Over Trials" 2nd edition.
Chapman & Hall, London.

Bing-Shun Wang, Xiao-Jin Wang, Li-Kun Gong (2009)
"The Construction of a Williams Design and Randomization
in Cross-Over Clinical Trials Using SAS"
J. of Statistical Software, Volume 29, Code Snippet 1
}
\author{
D. Labes
}
\note{
If ntmt>3 the design returned is randomly chosen from the possible designs.\cr\cr
A similar but simpler function was contained in the package \code{crossdes}.
But this package was removed from the CRAN repository.
}

\seealso{
\code{\link{sequences}}, \code{\link{RL4}}
}
\examples{
# Williams design for 4 treatments in ABC... notation
williams()
# The 6 sequences of the Williams design for 3 treatments
# same as sequences(design="3x6x3") except the ordering
williams(ntmt=3)
}
