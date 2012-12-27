\name{RL4}
\alias{RL4}
\title{
(Block) randomization of subjects, patients into sequence groups
}
\description{
This function creates a randomization list of subjects, patients 
into sequences of treatments.\cr
It is mainly intended for use in crossover studies but may be used also
for parallel group designs if for sequences c("T","R") is used.
}
\usage{
RL4(nsubj, seqs = c("TR", "RT"), blocksize, seed=runif(1,max=1E7), 
    randctrl=TRUE, pmethod=c("normal", "exact", "cc"), alpha=0.025)
}
\arguments{
  \item{nsubj}{
Number of subjects, patients ... to be randomized.\cr
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
  \item{pmethod}{
A character string describing the method for the p-value calculation of the runs
test.\cr\cr
\code{"normal"} denotes the normal approximation like the function(s) \code{runs.test()}
of the packages \code{tseries} or \code{lawstat}.\cr
\code{"exact"} chooses the calculation via exact distribution of the # of runs.\cr
\code{"cc"} chooses the continuity correction to the large sample approximation 
like in the statistical software SPSS.\cr\cr
Default is \code{pmethod="normal"}. 
}
  \item{alpha}{
Critical alpha value for re-creation of the randomlist based on the runs test.\cr
Recommended is \code{alpha=0.025}.
  }
}
\details{
As default the randomization is done as block randomization.\cr
If \code{blocksize} is given as a vector of numbers the sizes of the blocks 
are randomly chosen from that vector.\cr
If \code{blocksize=0} then the randomization takes place in the one block 
with all subjects.\cr\cr
The randomization is in the default settings controlled by a Wald-Wolfowitz
runs test of randomness, i.e. if this test shows significant deviation from randomness
(p-value<0.025 as default) the randomization list is recreated.\cr 
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
If this is not the case a warning is thrown that the design is not balanced.\cr\cr
The default for \code{pmethod} is the calculation via standard normal approximation.
This has shown the 'best' potential in rejecting the randomness for alternating
random lists.\cr
The randomness control does not work properly if more than 2 sequences are used.
In that case a test of repeated patterns of sequences in blocks of length=nseq is
done additionally to the runs test. Example (seqno) 1 2 3 1 2 3 ... 
}
\seealso{
\code{\link{print.rl4}, \link{runs.pvalue}}
}
\examples{
# block randomization of 12 subjects for a 2x2 cross-over using the defaults
# seqs = c("TR", "RT"), blocksize=4 and seed from Sys.time()
RL4(nsubj=12)

# block randomization of a full replicate design with randomly
# chosen blocksizes 2 or 4
rl <- RL4(nsubj=12, blocksize=c(2, 4), seqs=c("TRRT","RTTR"))
print(rl)

# randomization of 12 subjects for a 2x2 cross-over
# where the original random list don't pass the runs test
# watch the message.
RL4(nsubj=12, blocksize=0, seed=10)
#if you don't worry about some 'extreme' lists
RL4(nsubj=12, blocksize=0, seed=10, randctrl=FALSE)

# randomness control does'nt work that good in case of more 
# than 2 sequences: 3x3 latin square example
rl <- RL4(nsubj=12, seq <- c("ABC","BCA","CAB"), blocksize=3, seed=5125)
print(rl, sumry=TRUE)
# gives runs p.value=0.2502
}
