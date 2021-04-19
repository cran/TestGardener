
make.dataList <- function(U, key, optList, scrrng=NULL, titlestr=NULL,
                          nbin=nbinDefault(N), NumBasis=NumBasisDefault(N),
                          jitterwrd=TRUE, PcntMarkers=c( 5, 25, 50, 75, 95)) {
  
#' This function sets up the information required to analyze a set of data.
#' The information is stored in the struct object dataStr.
#' The information set up here is not meant to be modified by later code
#' in the analysis.
  
#  Last modified 18 March 2021 by Jim Ramsay

N <- nrow(U)
n <- ncol(U)

##  Set up UList

noption <- matrix(0,n,1)
for (item in 1:n) {
    noption[item] <- length(optList$optScr[[item]])# JL 2021-02-18
}
Wdim   <- sum(noption) # JL 2021-02-17

# Timings suggest U2List not worth the trouble.
# UList <- U2List(U, noption)

## compute sum scores for (both examinees and items

scrvec <- matrix(0,N,1)
itmvec <- matrix(0,n,1)
#Wdim   <- 0 # JL 2021-02-17
for (i in 1:n) {
  #Wdim <- Wdim + length(optList[[i]])
  for (j in 1:N) {
    scoreij   <- optList$optScr[[i]][U[j,i]]
    if (length(scoreij) > 0) {
      if (is.na(scoreij)) print(paste("is.na score:",j))
      if (is.null(scoreij)) print(paste("is.null score:",j))
    } else {
      print(paste("score of length 0:",j))
    }
    scrvec[j] <- scrvec[j] + scoreij
    itmvec[i] <- itmvec[i] + scoreij
  }
}

scrmin  <- min(scrvec)
scrmax  <- max(scrvec)
if (is.null(scrrng)) {
  scrrng <- c(scrmin,scrmax)
}
nfine   <- 101
scrfine <- seq(scrrng[1],scrrng[2],len=nfine)
denscdf <- seq(0,1,len=nfine)

thetaQnt <- seq(0,100,len=2*nbin+1)

##  jitter sum scores if required

if (jitterwrd) {
    scrjit <- scrvec + rnorm(N)*0.1
    scrjit[scrjit < scrmin] <- scrmin
    scrjit[scrjit > scrmax] <- scrmax
} else {
    scrjit <- scrvec
}

##  compute ranks for jittered sum scores

scrrnk <- matrix(0,N,1)
for (j in 1:N) {
    scrrnk[j] <- sum(scrjit <= scrjit[j])
}

percntrnk <- 100*scrrnk/N

##  Basis and bin setup for W function and theta estimation cycle

#  FdPar object for representing functions
#  The order of the B-splines is 5 because we need a 
#  smooth first derivative.

Wnorder <- 5  #  Order of the basis functions
# Set up the basis object
Wbasis  <- fda::create.bspline.basis(c(0,100), NumBasis, Wnorder) 
Wlambda <- 0   #  smoothing parameter
#  Compute the penalty matrix for the third derivative
Wnderiv <- 3  
Wpenmat <- fda::eval.penalty(Wbasis, Wnderiv)
#  Assemble this information into a fdPar object.
WfdPar  <- fdPar(Wbasis, Wnderiv, Wlambda, TRUE, Wpenmat)

##  Set up chartList to contain zero initial values for surprisal smoothing

#  This cell object is meant to be extracted and copied prior to the
#  analysis of the data so that initial values for the surprisal chart
#  may be changed over iterations.

chartList <- Wbinsmth.init(percntrnk, nbin, WfdPar, optList, U) 
  
##  Construct dataList object to define data Listucture

dataList <- list(U           = U, 
                 optList     = optList, 
                 key         = key,
                 chartList   = chartList,
                 WfdPar      = WfdPar, 
                 noption     = noption, 
                 nbin        = nbin, 
                 scrrng      = scrrng, 
                 scrfine     = scrfine,
                 scrvec      = scrvec, 
                 itmvec      = itmvec, 
                 percntrnk   = percntrnk, 
                 thetaQnt    = thetaQnt,
                 Wdim        = Wdim, 
                 PcntMarkers = PcntMarkers,
                 titlestr    = titlestr)

return(dataList)

}

#  ---------------------------------------------------------------

nbinDefault <- function(N) {
  if (N <= 500)              nbin <- floor(N/25)  
  if (N >  500 && N <= 2000) nbin <- floor(N/50)  
  if (N > 2000 && N <= 1e4)  nbin <- floor(N/100) 
  if (N >  1e4)              nbin <- 100 
  return(nbin)
}
  
#  ---------------------------------------------------------------

NumBasisDefault <- function(N) {
  if (N <= 500)               NumBasis <- 7                          
  if (N >  500 && N <= 2000)  NumBasis <- round(-14.7 + 8*log10(N))  
  if (N > 2000 && N <= 1e4)   NumBasis <- round(-14.7 + 8*log10(N))  
  if (N >  1e4)               NumBasis <-  24                        
  return(NumBasis)
}

#  ---------------------------------------------------------------

Wbinsmth.init <- function(percntrnk, nbin, WfdPar, optList, U) {
  
  # Last modified 22 February 2021 by Jim Ramsay
  
  #  This version of Wbinsmth() uses direct least squares smoothing of the
  #  surprisal values at bin centers to generate dependent variables for
  #  a linear model for the vectorized K by M-1 parameter matrix Bmat.
  
  nitem <- ncol(U)
  chartList <- vector("list", nitem)
  indfine  <- seq(0,100, len=101)
  thetaQnt <- seq(0,100, len=2*nbin+1)  
  bdry     <- thetaQnt[seq(1,2*nbin+1,by=2)]
  aves     <- thetaQnt[seq(2,2*nbin,  by=2)]  
  freq <- matrix(0,nbin,1)
  freq[1] <- sum(percntrnk < bdry[2])
  for (k in 2:nbin) {
    freq[k] <- sum(bdry[k-1] < percntrnk & percntrnk <= bdry[k])
  }
  meanfreq <- mean(freq)
  #  set up objects for required defining WfdPari for each item
  Wbasis    <- WfdPar$fd$basis
  Wnbasis   <- Wbasis$nbasis
  WLfd      <- WfdPar$Lfd
  Wlambda   <- WfdPar$lambda
  Westimate <- WfdPar$estimate
  Wpenmat   <- WfdPar$penmat
  WfdList   <- vector("list",nitem)
  for (item in 1:nitem) {
    # print(item)
    Mi    <- length(optList$optScr[[item]])
    logMi <- log(Mi)
    Uveci <- as.numeric(U[,item])
    Pbin  <- matrix(0,nbin,Mi)  #  probabilities
    Wbin  <- matrix(0,nbin,Mi)  #  transformation of probability
    indpts <- c(1:nbin)
    for (k in 1:nbin) {
      #  index of percntrnk values within this bin
      indk   <- percntrnk >= bdry[k] & percntrnk <= bdry[k+1]
      if (sum(indk) > 0) {
        Uvecik <- Uveci[indk]
        nk     <- sum(indk)
        for (m in 1:Mi) {
          Pbin[k,m] <- sum(Uvecik == m)/nk
          if (Pbin[k,m] == 0) Pbin[k,m] <- NA
        }
        Wbin[k,] <- -log(Pbin[k,])/logMi
      } else {
        Pbin[k,] <- NA
      }
    } # end of bin loop
    #  unrestricted smoothing for (the multi-option case
    maxWbin <- 0
    for (m in 1:Mi) {
      Wmis.na <- is.na(Pbin[,m])
      indm <- (1:nbin)[!Wmis.na]
      if (length(indm) > 0) maxWbin <- max(c(maxWbin,max(Wbin[indm,m])))
    }
    SurprisalMax <- min(c(-log(1/(meanfreq*2))/logMi, maxWbin))
    for (m in 1:Mi) {
      Wmis.na <- is.na(Pbin[,m])
      indm <- (1:nbin)[!Wmis.na]
      # if (m < Mi) {
      #   Wbin[Wmis.na,m] <- SurprisalMax
      # } else {
        indmlen <- length(indm)
        if (indmlen > 3) {
          WY <- Wbin[indm,m];
          WX <- cbind(rep(1,indmlen), aves[indm])
          BX <- lsfit(aves[indm], WY)$coefficients
          Wbin[indm,m] <- WX %*% BX
          Wbin[Wmis.na,m] <- SurprisalMax
        } else {
          Wbin[Wmis.na,m] <- SurprisalMax
        # }
      }
    }
    #  generate a map into M-vectors with zero row sums
    if (Mi == 2) {
      root2 <- sqrt(2)
      Zmati <- matrix(1/c(root2,-root2),2,1)
    } else {
      Zmati <- zerobasis(Mi)
    }
    
    #  apply conventional surprisal smoothing 
    Wfdi     <- fd(matrix(0,Wnbasis,Mi-1),Wbasis)
    WfdPari  <- fdPar(Wfdi, WLfd, Wlambda, Westimate, Wpenmat)
    Sfdi     <- fda::smooth.basis(aves, Wbin, WfdPari)$fd
    #  compute spline basis functions at bin centres
    Phimati <- fda::eval.basis(aves, WfdPar$fd$basis)
    #  evaluate smooth at bin centres
    Smathati <- fda::eval.fd(aves, Sfdi)
    #  map this into zero-row-sum space
    Smatctri <- Smathati %*% Zmati
    #  regress the centred data on the negative of basis values
    Result   <- lsfit(-Phimati, Smatctri, intercept=FALSE)
    SSE0 <- sum(Smatctri^2)
    SSE1 <- sum(Result$residuals^2)
    RSQR <- (SSE0-SSE1)/SSE0
    Bmathati <- Result$coeff
    chartList[[item]] <- matrix(Bmathati,WfdPar$fd$basis$nbasis,Mi-1)
  }
  
  return(chartList)
  
}

    
    
