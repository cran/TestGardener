TG_analysis <- function(chcemat, scoreList, noption, sumscr_rng=NULL, 
                        titlestr=NULL, itemlabvec=NULL, optlabList=NULL,
                        nbin=nbinDefault(N), NumBasis=7, SfdPar=NULL,
                        jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
                        ncycle=10, cyclechoice=ncycle, 
                        itdisp=FALSE, verbose=FALSE) {
  
  #. Arguments:
  #. chcemat     ... An N by n matrix.  Column i must contain the integers 
  #.                 from 1 to M_i, where M_i is the number of options
  #.                 for item i.  If missing or illegitimate responses exist
  #.                 for item i,  the column must also contain an integer
  #.                 greater than M_i that is used to identify such responoses. 
  #.                 Alternatively, the column use NA for this purpose.
  #.                 Because missing and illegible responses are normally
  #.                 rare, they are given a different and simpler estimation
  #.                 procedure for their surprisal values.
  #.                 chcemat is mandatory.
  #. scoreList   ... Either a list of length n, each containing a vector of 
  #.                 length M_i that assigns numeric weights to the options
  #.                 for that item.  
  #.                 In the special case of multiple choice items where the 
  #.                 correct option has weight 1 and all others weight 0, 
  #.                 a single integer can identify the correct answer.
  #.                 If all the items are of the multiple 
  #.                 type, scoreList may be a numeric vector of length n
  #.                 containing the right answer indices.  List object
  #.                 scoreList is mandatory because these weights define the
  #.                 person scores for the surprisal curve estimation process.
  #. noption     ... A numeric vector containing the number of choices for each
  #.                 item.  These should not count missing or illegal choices.
  #.                 Although this object might seem redundant, it is needed
  #.                 for checking the consistencies among other objects and
  #.                 as an aid for detecting missing and illegal choices.
  #. sumscr_rng  ... A vector of length 2 indicating the initial and final
  #.                 sum score values.  Default is NULL the whole sum score
  #.                 is used.
  #. titlestr    ... A title string for the data and their analyses.
  #.                 Default is NULL.
  #. itemlabvec  ... A character value containing labels for the items.
  #.                 Default is NULL and item position numbers are used.
  #. optlabList. ... A list vector of length n, each element i of which is a
  #.                 character vector of length M_i.
  #.                 Default is NULL, and option numbers are used.
  #. nbin        ... The number of bins containing proportions of choices.
  #. NumBasis    ... The number of spline basis functions to use for 
  #.                 surprisal values.  Defaults to 7.
  #. SfdPar      ... A functional parameter object to be used for representing
  #.                 surprisal curves.  May be NULL in which case the spline  
  #                  basis has NumBasis basis functions and is of order 5.
  #. jitterwrd   ... A logical object indicating whether a small jittering
  #.                 perturbation should be used to break up ties.  
  #                  Defaults to TRUE.
  #. PcntMarkers ... A vector of percentages inside of [0,100] that appear
  #                  in plots.  Defaults to c(5, 25, 50, 75, 95).
  #. verbose     ... Extra displays are provided.  Defaults to FALSE.
  #. ncycle      ... The number of cycles in the analysis.  Defaults to 10.
  #. choicecycle ... A number within 1 to 10 indicating which cycle will be 
  #.                 used to represent the TestGardener results.  
  #.                 Defaults to ncycle.

    #  Last modified 3 November 2023 by Jim Ramsay
  
  N <- nrow(chcemat)
  n <- ncol(chcemat)
  
  dataList <- make_dataList(chcemat, scoreList, noption, sumscr_rng=NULL, 
                            titlestr, itemlabvec, optlabList,
                            nbin, NumBasis, SfdPar,
                            jitterwrd, PcntMarkers)
  
  #  ----------------------------------------------------------------------------
  #  compute the initial option surprisal curves using the 
  #  percentage ranks as initial estimates of index
  #  ----------------------------------------------------------------------------
  
  index    <- dataList$percntrnk
  indexQnt <- dataList$indexQnt
  
  #  ----------------------------------------------------------------------------
  #                      Proceed through the cycles
  #  ----------------------------------------------------------------------------
  
  parmListvec <- Analyze(index, indexQnt, dataList,   
                         ncycle, itdisp, verbose) 
  
  # print("Analyze complete")
  #  ----------------------------------------------------------------------------
  #              Plot the average H value, meanHsave, over cycles
  #  ----------------------------------------------------------------------------
  
  # print(ncycle)
  HALsave <- matrix(0,ncycle,2)
  for (icycle in 1:ncycle) {
    HALsave[icycle,1] <- parmListvec[[icycle]]$meanF
    HALsave[icycle,2] <- parmListvec[[icycle]]$infoSurp
  }
  
  par(mfrow=c(2,1))
  
  plot(1:ncycle, HALsave[,1], type="b", lwd=2, 
       xlab="Cycle Number",ylab="Mean H")
  plot(1:ncycle, HALsave[,2], type="b", lwd=2, 
       xlab="Cycle Number", ylab="Arc Length")
  
  # print("cycle plotting complete")
  
  #  ----------------------------------------------------------------------------
  #             Select cycle for results in parmList cnd define it
  #  ----------------------------------------------------------------------------
  
  parmList  <- parmListvec[[cyclechoice]]
  
  #  set up estimated objects after iterations
  
  SfdList   <- parmList$SfdList
  Qvec      <- parmList$Qvec
  binctr    <- parmList$binctr
  index     <- parmList$index
  infoSurp  <- parmList$infoSurp
  
  # print("transfer from parmList complete")
  
  #  ----------------------------------------------------------------------------
  #   Compute the arc length or information measure and objects it needs for 
  #   plotting results as a function of arc length
  #  ----------------------------------------------------------------------------
  
  # print("entering index2info")
  infoList <- index2info(index, Qvec, SfdList, binctr)
  # print("index2info complete")

  return(list(dataList=dataList, parmList=parmList, infoList=infoList))

}

#  ---------------------------------------------------------------

nbinDefault <- function(N) {
  if (N <= 500)              nbin <- floor(N/25)  
  if (N >  500 && N <= 2000) nbin <- floor(N/50)  
  if (N > 2000 && N <= 1e4)  nbin <- floor(N/100) 
  if (N >  1e4)              nbin <- 100 
  return(nbin)
}
