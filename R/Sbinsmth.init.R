#  ---------------------------------------------------------------

Sbinsmth.init <- function(percntrnk, nbin, Sbasis, grbgvec, noption, chcemat) {
  
  # Last modified 24 March 2024 by Jim Ramsay
  
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


