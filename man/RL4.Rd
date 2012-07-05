\name{RL4}
\alias{RL4}
\title{
(Block) randomisation of subjects, patients into sequence groups
}
\description{
This function creates a randomisation list of subjects, patients 
into sequences of treatments.\cr
It is mainly intended for use in crossover studies but may be used also
for parallel group designs if for sequences c("T","R") is used.
}
\usage{
RL4(nsubj, seqs = c("TR", "RT"), blocksize, seed = Sys.time(), 
    randctrl=TRUE, alpha=0.025)
}
\arguments{
  \item{nsubj}{
Number of subjects, patients ... to be randomised.\cr
Or a vector of the subject numbers (f.i. 1001:1012)
  }
  \item{seqs}{
Character representation of the sequences.\cr
In BE studies for a 2x2 cross-over usually something like \code{c("TR", "RT")}.\cr
If you prefer the ABC notation write down here f.i. \code{c("AB", "BA")}.
  }
  \item{blocksize}{
Size of the blocks to randomize within.\cr
If blocksize is missing it defaults to 2*number of sequences.\cr\cr
\code{blocksize} may be a vector of numbers. \cr
In that case the sizes of the blocks are randomly chosen from that vector.\cr\cr
If \code{blocksize=0} then the randomization takes place in the one block 
with all subjects.
  }
  \item{seed}{
An integer for the initialisation of the random number generator.
}
  \item{randctrl}{
Shall the creation of the randomlist controlled by a runs test of randomness?\cr
Defaults to \code{TRUE}.\cr
  }
  \item{alpha}{
Critical alpha value for re-creation of the randomlist based on the runs test.\cr
Recommended is \code{alpha=0.025}.
  }
}
\details{
As default the randomisation is done as block randomisation.\cr
If \code{blocksize} is given as a vector of numbers the sizes of the blocks 
are randomly chosen from that vector.\cr
If \code{blocksize=0} then the randomization takes place in the one block 
with all subjects.\cr\cr
The randomisation is in the default settings controlled by a Wald-Wolfowitz
runs test of randomness, i.e. if this test shows significant deviation from randomness
(p-value<0.025 as default) the randomisation list is recreated.\cr 
This behavior may be switched off by setting \code{randctrl=FALSE} if you don't
see any needfulness for such a control.\cr
The reason for such a control was originally to avoid randomlists with only
'alternating' sequences f.i. "TR" "RT" "TR" "RT" ...\cr
See \url{http://forum.bebac.at/mix_entry.php?id=8745}\cr\cr
In its original form the runs test is only defined for dichotomous outcomes 
(i.e if 2 sequences are used).\cr
If more than 2 sequences are used the runs test is modified by dichotomizing
the sample of sequence numbers into cases < median and > median.\cr
See package \code{lawstat} or \url{http://de.wikipedia.org/wiki/Run-Test} for this.
}
\value{
Returns a list of class "rl4".
This list has the components 
  \item{rl }{A data.frame with the random list.}
  \item{seed }{The used seed.}
  \item{blocksize }{The used blocksize.}
  \item{ninseqs }{A named vector containing the number of subjects in the sequence groups.}
  \item{runs.pvalue }{The p-value of a runs test of randomness.}
  \item{date }{Date and time of creation.}
\cr
The class rl4 has the S3 method \code{\link{print.rl4}}.
}

\author{
D. Labes\cr
Part of the code for using the runs test for randomness 
according to a suggestion by Helmut Schuetz.
}
\note{
The blocksize(s) should be a multiple of the used # of sequences. 
If this is not the case it (they) will be adjusted accordingly.\cr\cr
The number of subjects should be a multiple of the used # of sequences. 
If this is not the case a warning is thrown that the design is not balanced.
}
\seealso{
\code{\link{print.rl4}, \link{runs.pvalue}}
}
\examples{
# block randomisation of 12 subjects for a 2x2 cross-over using the defaults
# seqs = c("TR", "RT"), blocksize=4 and seed from Sys.time()
RL4(nsubj=12)

# block randomisation of a full replicate design with randomly
# chosen blocksizes 2 or 4
rl <- RL4(nsubj=12, blocksize=c(2, 4), seqs=c("TRRT","RTTR"))
print(rl)

# randomisation of 12 subjects for a 2x2 cross-over
# where the original random list don't pass the runs test
# watch the message.
RL4(nsubj=12, blocksize=0, seed=10)
#if you don't worry about some 'extreme' lists
RL4(nsubj=12, blocksize=0, seed=10, randctrl=FALSE)
}
