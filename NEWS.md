# Version 0.3-6

released 2023-08-19

* Man page of package removed because version and date became outdated
  (requested by Kurt Hornik).

# Version 0.3-5

released 2019-08-24

* Bug fix in function sequences() wrt design="parallel".
THX to Yung-jin Lee

# Version 0.3-3

released 2017-03-22

* Maintenance release according to a request of CRAN to avoid warnings /
errors in control statements with condition with length > 1

# Version 0.3-2

released 2015-06-30

* Maintenance release with NAMESPACE adapted to get the necessary functions from package stats 
(necessary for R devel)

# Version 0.3-1

released 2012-12-27

* Function williams() introduced to obtain the sequences of a Williams design

# version 0.2-1

released 2012-09-20

* Function sequences introduced to obtain the sequences of commonly used designs in BE studies

# version 0.1-3

released 2012-08-10

* Exact version of Wald-Wolfowitz runs test implemented. See ?runs.pvalue.
* Random list now contains also a sequence number. May be better suited 
to detect patterns in the random list if more than 2 sequences are used. 
* 'Randomness' control in function RL4() in case of more than 2 sequences now tests the occurrence of repeated patterns of sequences hard coded.
Example seqno= 1 2 3 1 2 3 ... will be recognized and avoided.

# version 0.1-2

released 2012-07-05

* First release via CRAN.
