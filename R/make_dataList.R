
make_dataList <- function(chcemat, scoreList, noption, sumscr_rng=NULL, 
                          titlestr=NULL, itemlabvec=NULL, optlabList=NULL,
                          nbin=nbinDefault(N), NumBasis=7, 
                          jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95),
                          verbose=FALSE) {
  
  #' This function sets up the information required to analyze a set of data.
  #' The information is stored in the struct object dataStr.
  #' The information set up here is not meant to be modified by later code
  #' in the analysis.
  
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
  #.                 Object chcemat is mandatory.
  #. scoreList:      Either a list of length n, each containing a vector of 
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
  #.                 Object scoreList is mandatory.
  #. noption     ... A numeric vector containing the number of choices for each
  #.                 item.  These should not count missing or illegal choices.
  #.                 Although this object might seem redundant, it is needed
  #.                 for checking the consistencies among other objects and
  #.                 as an aid for detecting missing and illegal choices.
  #.                 Object noptions is mandatory.
  #. sumscr_rng  ...  A vector of length 2 indicating the initial and final
  #.                  sum score values.  Default is NULL the whole sum score
  #.                  is used..
  #. titlestr    ...  A title string for the data and their analyses.
  #.                  Default is NULL.
  #. itemlabvec  ...  A character value containing labels for the items.
  #.                  Default is NULL and item position numbers are used.
  #. optlabList. ...  A list vector of length n, each element i of which is a
  #.                  character vector of length M_i.
  #.                  Default is NULL, and option numbers are used.
  #. nbin        ...  The number of bins containing proportions of choices.
  #. NumBasis    ...  The number of spline basis functions to use for 
  #.                  surprisal values.  Defaults to 7.
  #. jitterwrd   ...  A logical object indicating whether a small jittering
  #.                  perturbation should be used to break up ties.  
  #                   Defaults to TRUE.
  #. PcntMarkers ...  A vector of percentages inside of [0,100] that appear
  #                   in plots.  Defaults to c(5, 25, 50, 75, 95).
  #. verbose     ...  Extra displays are provided.  Defaults to FALSE.
  
  #  Last modified 16 September 2024 by Jim Ramsay
  
  # Default parameter values
  
  # sumscr_rng=NULL
  # titlestr=NULL 
  # itemlabvec=NULL 
  # optlabList=NULL
  # nbin=nbinDefault(N) 
  # NumBasis=7
  # jitterwrd=TRUE 
  # PcntMarkers=c( 5, 25, 50, 75, 95)
  # verbose=TRUE
  
  #. Dimensions of index matrix
  
  # print("inside make_dataList")
  
  N <- nrow(chcemat)
  n <- ncol(chcemat)
  
  #. --------------------------------------------------------------------------
  #.                          Check Arguments
  #. -------------------------------------------------------------------------- 
  
  #  Check index matrix chcemat (short for "choice matrix")
  
  if (is.null(chcemat)) stop("Index matrix chcemat is NULL")
  if (min(chcemat) < 1) stop("Non-positive index values encountered.")
  # if (any(!is.integer(chcemat))) stop("Non-integer values found in chcemat.")
  if (any(is.na(chcemat))) stop("One or more values in chcemat are NA")
  if (N < 100) warning(paste("Number of rows < 100,", 
                  "A TestGardener analysis is not likely to succeed."))
  
  #.    Check vector noption, containingthe number of options in each item
  #.    Garbage options indicating missing of illegal values are not added.

  if (is.null(noption))          stop("Vector noption is NULL.")
  if (length(noption) != n)      stop("Argument noption is not of length n.")
  if (any(is.na(noption)))       stop("A value in noption is NA.")
  if (min(noption) < 2)          stop("Values less than 2 in noption.")
  # if (any(!is.integer(noption))) stop("Non-integer values in noption.")
  
  # compute dimension of ambient space
  
  Sdim  <- sum(noption) 
  
  #  Check if scoreList is.numeric and the test is all multiple choice items
  key <- NULL
  if (is.numeric(scoreList)) {
    #. scoreList is in key format, numeric vector
    print(paste("Argument scoreList is numeric,",
               "and chcemat has multiple choice data."))
    key <- scoreList
    if (min(key) < 1) stop(paste("Zero data values encountered in key.",
                                 "Are they score values?"))
    if (max(abs(as.integer(key)-key)) > 0) 
      stop("Non-integer values found in key.")
    if (length(key) != n && length(key) > 0) 
      stop("length of key is neither n nor 0") 
    for (i in 1:n) {
      keyerror <- FALSE
      if (key[i] > noption[i]) {
        print(paste("key value larger than number of options for item",
                    i,"."))
        keyerror <- TRUE
      }
    }
    if (keyerror) stop("")
    #  Set up key as scoreList,  a numbered list.
    #. each list element is a vector with values 0 or 1
    scoreList <- vector("list",n)
    for (i in 1:n) {
      scoreveci  <- rep(0,noption[i])
      scoreveci[key[i]] <- 1
      scoreList[[i]] <- scoreveci
    }
  }
  
  #  check itemlabvec
  
  if (!is.null(itemlabvec)) {
    if (length(itemlabvec) != n ) stop("Number of item labels is not n.")
    if (!is.character(itemlabvec)) itemlabvec <- as.character(itemlabvec)
  }
  
  #  check optlabList
  
  if (!is.null(optlabList)) {
    if (length(optlabList) != n ) 
      stop("Number of option label sets is not n.")
    for (i in 1:n) {
      for (m in 1:noption[i]) 
        if (!is.character(optlabList[[i]])) 
          optlabList[[i]] <- as.character(optlabList[[i]])
    }
  }
  
  #  check sumscr_rng
  
  if (!is.null(sumscr_rng) & !is.numeric(sumscr_rng))
    stop("Argument sumscr_rng is neither NULL nor numeric.")
  if (is.numeric(sumscr_rng) & length(sumscr_rng) != 2) 
    stop("Argument sumscr_rng is not of length 2.")
  
  #. --------------------------------------------------------------------------
  #  Construct logical vector grbgvec indicating items with missing or 
  #  illegal choices and therefore needing an extra garbage option
  #. If so, update noption by one.
  #. --------------------------------------------------------------------------
  
  grbgvec <- rep(FALSE, n)
  for (item in 1:n) {
    grbgvecind <- chcemat[,item] > noption[item]
    if (sum(grbgvecind) > 0) {
      #  indices greater than noption[item] encountered
      #  add an option the these labels, set grbgvec value to TRUE
      #  change indices to updated noption[item]
      noption[item]            <- noption[item] + 1
      grbgvec[item]            <- TRUE
      chcemat[grbgvecind,item] <- noption[item]
      if (!is.null(optlabList))
          optlabList[[item]] <- c(optlabList[[item]],0)
    }
  }
  
  #. --------------------------------------------------------------------------
  #         Compute sum scores for both examinees and items
  #. --------------------------------------------------------------------------
  
  scrvec <- matrix(0,N,1)
  itmvec <- matrix(0,n,1)
  for (item in 1:n) {
    for (j in 1:N) {
      scorevec <- scoreList[[item]]
      scoreij  <- scorevec[chcemat[j,item]]
      if (!is.na(scoreij)) { 
        scrvec[j] <- scrvec[j] + scoreij
        itmvec[item] <- itmvec[item] + scoreij
      }
    }
  }
  
  #  set up a fine mesh of score index values over sum score range
  
  scrmin  <- min(scrvec)
  scrmax  <- max(scrvec)
  if (is.null(sumscr_rng)) sumscr_rng <- c(scrmin,scrmax)
  nfine   <- 101
  scrfine <- seq(sumscr_rng[1],sumscr_rng[2],len=nfine)
  
  #  initial set of bin centres and boundaries
  
  indexQnt <- seq(0,100,len=2*nbin+1)
  
  #. --------------------------------------------------------------------------
  #  jitter sum scores to break up ties in integer valued sum scores
  #  if jitterwrd is TRUE, the jitter values are set up here
  #  if jitterwrd is FALSE, sum scores are not jittered
  #  if jitterwrd is a numeric vector of length N containing jittered values
  #    these are used.  This allows multiple data analyses without
  #    change of jittered sum score values.
  #  if jitterwrd is not logical and is not a numeric vector of length N
  #    an error is declared.
  #. --------------------------------------------------------------------------
  
  if (is.logical(jitterwrd)) {
    if (jitterwrd) {
      jitter <- rnorm(N)*0.1
    } else {
      jitter <- rep(0,N)
    }
    scrjit <- scrvec + jitter
    scrjit[scrjit < scrmin] <- scrmin
    scrjit[scrjit > scrmax] <- scrmax
  } else {
    if (is.numeric(jitterwrd) && length(jitterwrd) == N) {
      scrjit <- jitterwrd
      scrjit[scrjit < scrmin] <- scrmin
      scrjit[scrjit > scrmax] <- scrmax
    } else {
      stop(paste("jitterwrd is not logicial",
                 "and is not a numeric vector of length N"))
    }
  }
  
  #  compute ranks for jittered sum scores
  
  scrrnk <- matrix(0,N,1)
  for (j in 1:N) scrrnk[j] <- sum(scrjit <= scrjit[j])
  percntrnk <- 100*scrrnk/N
  
  #. --------------------------------------------------------------------------
  ##     Set up spline basis and bins for initial surprisal curves.
  #. --------------------------------------------------------------------------
  
  #  number of basis functions.  If NULL, this is assigned according to size of N
  
  if (round(NumBasis) == NumBasis) {
    #  NumBasis is an integer
    if (NumBasis < 2) 
      stop("There must be at least two spline basis functions.")
    if (NumBasis > 7) 
      warning("More than 7 basis functions may cause instability 
                in optimization.")
    Snbasis <- NumBasis
    Snorder <- min(Snbasis, 5)
    Sbasis  <- fda::create.bspline.basis(c(0,100), Snbasis, Snorder)
  } else {
    stop("Number of basis functions is not an integer.")
  }
  
  #. --------------------------------------------------------------------------
  #     Sbinsmth.init computes initial approximations to surprisal curves
  #.    using functional data smoothing.
  #. --------------------------------------------------------------------------
  
  # print("entering Sbinsmth.init")
  SfdList <- Sbinsmth.init(percntrnk, nbin, Sbasis, grbgvec, noption, chcemat) 
  # print("left Sbinsmth.init")
  
  ##  Construct dataList object to define data Listucture
  
  dataList <- list(chcemat     = chcemat, 
                   scoreList   = scoreList,
                   noption     = noption, 
                   titlestr    = titlestr,
                   itemlabvec  = itemlabvec,
                   optlabList  = optlabList,
                   SfdList     = SfdList,
                   key         = key,
                   grbgvec     = grbgvec,
                   nbin        = nbin, 
                   sumscr_rng  = sumscr_rng, 
                   scrfine     = scrfine,
                   scrvec      = scrvec,
                   scrjit      = scrjit,
                   itmvec      = itmvec, 
                   percntrnk   = percntrnk, 
                   indexQnt    = indexQnt,
                   Sdim        = Sdim, 
                   PcntMarkers = PcntMarkers,
                   N           = N,
                   n           = n,
                   NumBasis    = NumBasis)
  
  return(dataList)
  
}

#  ---------------------------------------------------------------

Sbinsmth.init <- function(percntrnk, nbin, Sbasis, grbgvec, noption, chcemat) {
  
  # Last modified 26 November 2023 by Jim Ramsay
  
  #  This version of Sbinsmth() uses direct least squares smoothing of the
  #  surprisal values at bin centers to generate dependent variables for
  #  a model for the vectorized K by M-1 parameter matrix Bmat.
  
  nitem     <- ncol(chcemat)
  chartList <- vector("list", nitem)
  indfine   <- seq(0,100, len=101)
  indexQnt  <- seq(0,100, len=2*nbin+1)  
  bdry      <- indexQnt[seq(1,2*nbin+1,by=2)]
  binctr    <- indexQnt[seq(2,2*nbin,  by=2)]  
  freq      <- matrix(0,nbin,1)
  freq[1]   <- sum(percntrnk < bdry[2])
  for (k in 2:nbin) {
    freq[k] <- sum(bdry[k-1] < percntrnk & percntrnk <= bdry[k])
  }
  meanfreq <- mean(freq)
  SfdList  <- vector("list", nitem)
  Sbasis   <- Sbasis
  Snbasis  <- Sbasis$nbasis
  
  #. --------------------------------------------------------------------------
  #. loop through items to compute their surprisal curve approximations
  #. --------------------------------------------------------------------------
  
  # print("entering loop")
  
  for (item in 1:nitem) {
    # print(paste("item =",item))
    Mi    <- noption[item]
    logMi <- log(Mi)
    chcematveci <- as.numeric(chcemat[,item])
    Pbin  <- matrix(0,nbin,Mi)  #  probabilities
    Sbin  <- matrix(0,nbin,Mi)  #  transformation of probability
    # print("k loop")
    for (k in 1:nbin) {
      #  index of percntrnk values within this bin
      indk   <- percntrnk >= bdry[k] & percntrnk <= bdry[k+1]
      if (sum(indk) > 0) {
        chcematvecik <- chcematveci[indk]
        nk     <- sum(indk)
        for (m in 1:Mi) {
          Pbin[k,m] <- sum(chcematvecik == m)/nk
          if (Pbin[k,m] == 0) Pbin[k,m] <- NA
        }
        Sbin[k,] <- -log(Pbin[k,])/logMi
      } else {
        Pbin[k,] <- NA
      }
    } # end of bin loop
    
    #  Smooth the binned S values
    
    #  Set up SurprisalMax to replace NA's
    
    # print("SurprisalMax")
    
    maxSbin <- 0
    for (m in 1:Mi) {
      Smis.na <- is.na(Pbin[,m])
      indm <- (1:nbin)[!Smis.na]
      if (length(indm) > 0) maxSbin <- max(c(maxSbin,max(Sbin[indm,m])))
    }
    SurprisalMax <- min(c(-log(1/(meanfreq*2))/logMi, maxSbin))
    
    #  process NA values in Sbin associated with zero probabilities
    
    # print("m loop")
    
    for (m in 1:Mi) {
      Smis.na <- is.na(Pbin[,m])
      if (!grbgvec[item] || (grbgvec[item] && m != Mi)) {
        Sbin[Smis.na,m] <- SurprisalMax
      }  else {
        #  garbage choices: compute sparse numeric values into 
        #  linear approximations and NA values to SurprisalMax
        indm    <- (1:nbin)[!Smis.na]
        indmlen <- length(indm)
        nonindm <- (1:nbin)[Smis.na]
        if (indmlen > 3) {
          SY <- Sbin[indm,m];
          SX <- cbind(rep(1,indmlen), binctr[indm])
          BX <- lsfit(binctr[indm], SY)$coefficients
          Sbin[indm,m]    <- SX %*% BX
          Sbin[nonindm,m] <- SurprisalMax
        } else {
          Sbin[nonindm,m] <- SurprisalMax
        }
      }
    }
    
    #  generate a map into M-vectors with zero row sums
    
    if (Mi == 2) {
      root2 <- sqrt(2)
      Zmati <- matrix(1/c(root2,-root2),2,1)
    } else {
      Zmati <- fda::zerobasis(Mi)
    }
    
    # print("smooth")
    # print(class(Sbin))
    #  apply conventional smoothing of surprisal values
    # print(class(Sbasis))
    Sfdi     <- fda::smooth.basis(binctr, Sbin, Sbasis)$fd
    # print(class(Sfdi))
    #  compute spline basis functions at bin centres
    Phimati  <- fda::eval.basis(binctr, Sbasis)
    #  evaluate smooth at bin centres
    Smathati <- fda::eval.fd(binctr, Sfdi)
    #  map this into zero-row-sum space
    Smatctri <- Smathati %*% Zmati
    #  regress the centred data on the negative of basis values
    Result <- lsfit(-Phimati, Smatctri, intercept=FALSE)
    Bmati  <- Result$coefficient
    Sfdi   <- fda::fd(Bmati, Sbasis)
    
    #  store objects in SListi
    
    # print("returned objects list")
    SListi <- list(
      Sfd        = Sfdi,       #  functional data object for (options
      M          = Mi,         #  the number of options
      Pbin       = Pbin,       # proportions at each bin
      Sbin       = Sbin,       # negative surprisals at each bin
      Zmat       = Zmati,
      Pmatfine   = NULL,   
      Smatfine   = NULL,   
      DSmatfine  = NULL,  
      D2Smatfine = NULL
    )
    SfdList[[item]] <- SListi
  }

  return(SfdList)
  
}

#  ---------------------------------------------------------------

nbinDefault <- function(N) {
  if (N <= 500)              nbin <- floor(N/25)  
  if (N >  500 && N <= 2000) nbin <- floor(N/50)  
  if (N > 2000 && N <= 1e4)  nbin <- floor(N/100) 
  if (N >  1e4)              nbin <- 100 
  return(nbin)
}


  



    
    
