###############################################################################
# function RL4 with randomness control via runs test
# 
# Author: dlabes
###############################################################################
# does not require(tseries) any longer
# function runs.pvalue implemented

# function for (block) randomization
RL4 <- function(nsubj, seqs=c("TR","RT"), blocksize, seed=Sys.time(),
                randctrl=TRUE, alpha=0.025)
{
  if (length(nsubj)>1) {
    subj <- nsubj
    nsubj <- length(subj)
  } else {
    subj <- 1:nsubj
  }
  seed <- as.numeric(seed)
  if(is.na(seed)) seed <- 0
  set.seed(seed)

  nseq <- length(seqs)
  
  if(missing(blocksize)) blocksize <- 2*nseq
    
  if(blocksize[1]>nsubj) {
    blocksize <- 0
    warning("Blocksize > # of subjects!", 
        " Blocksize adapted to ", nsubj,".", call. = FALSE )
  }
  
  
  if (blocksize==0) {
    blocksize <- nsubj
  } else {
    # silent adaption or warning?
    blocksize[blocksize<nseq] <- nseq 
    # is blocksize a multiple of sequences?
    # if not, a highly unbalanced design may be the result
    # therefore we adapt the blocksize to a multiple if # of seqs
    if (any(nseq*(blocksize%/%nseq) != blocksize)) {
      blocksize <- nseq*(blocksize%/%nseq)
      # may be that due to truncation some elements are <nseq
      blocksize[blocksize<nseq] <- nseq # vector form
      warning("Blocksize is not a multiple of sequences!", 
              " Blocksize adapted to ", blocksize,".", call. = FALSE )
      blocksize <- unique(blocksize)
    }
  }  
  # call the helper function
  rlv <- rlv(nsubj, nseq, blocksize)
  rl  <- rlv$rl
  # runs test of randomnes is only possible if 2 sequences
  # if more than 2 sequences we use the dichotomization by the median
  # see wikipedia entry or package lawstat
  runs.p <- runs.pvalue(rl)
  # randomness control
  if (randctrl){
    while(runs.p < alpha){
      msg <- paste("runs.p=", 
                   format(runs.p, digits=4),". Recreating randomlist.")
      message(msg)
      rlv <- rlv(nsubj, nseq, blocksize)
      rl  <- rlv$rl
      runs.p <- runs.pvalue(rl)
    }
  }
  bsv <- rlv$bsv
  # create the randomlist with character representation of seqs
  rlc <- seqs[rl]
  # check if design is balanced
  ns <- table(rlc)
  nsequal <- (ns - ns[1]) == 0
  if (!all(nsequal)){
    msg  <- paste(" ", names(ns))
    msg2 <- paste(" ", ns)
    warning("Unbalanced design!", " # of subj. in sequences", 
             msg, ":", msg2, call. = FALSE)
  }
  # the random list itself
  rl <- data.frame(subject=subj, sequence=rlc, stringsAsFactors=FALSE)
  # number of subjects in groups
  ns    <- t(as.matrix(ns))
  nsv   <- as.vector(ns)
  names(nsv) <- colnames(ns)
  # blocksize(s)
  if (length(unique(bsv))==1) bsv <- bsv[1]
  rlret <- list(rl=rl, seed=seed, blocksize=bsv, 
                ninseqs=nsv, runs.pvalue=runs.p, date=Sys.time())
  class(rlret) <- "rl4"
  return(rlret)
}

# internal function for (block) randomisation
# returns a list with the components
#   rl = random list (sequences numeric 1:nseq)
#   bsv = blocksize vector (actual)
rlv <- function(nsubj, nseq, blocksize)
{
  # we are working with the numeric coding of sequences
  seqsn <- 1:nseq
  n <- 0
  rl  <- vector(mode="numeric") # numeric random sequence 
  bsv <- vector(mode="numeric") # blocksize vector
  # loop over blocks
  while (n<nsubj){
    # choose a blocksize by random
    # if blocksize has only one element the sample function
    # chooses integers up to blocksize
    bs <- ifelse(length(blocksize)>1, sample(blocksize,1), blocksize[1])
    # last block may be smaller
    bs <- ifelse((nsubj-n) < bs, nsubj-n, bs)
    # random list of a block
    rpp <- ifelse(nseq*bs%/%nseq != bs, bs%/%nseq+1, bs%/%nseq)
    rlb <- sample(rep(seqsn, rpp), bs)
    # debug prints
    #cat("bs:");print(bs)
    #cat("rlb");print(rlb)
    rl <- c(rl, rlb)
    bsv <- c(bsv, bs)
    n <- n + bs
  }
  rl <- rl[1:nsubj]
  return(list(rl=rl, bsv=bsv))
}

