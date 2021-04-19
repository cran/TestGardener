smooth.surp <- function(argvals, Wbin, Bmat0, WfdPar, wtvec=NULL, conv=1e-4,
                       iterlim=50, dbglev=0) {
  #  Smooths the relationship of Y to ARGVALS using weights in WTVEC by fitting 
  #     surprisal functions to a set of surprisal transforms of choice 
  #     probabilities, where the surprisal transformation of each probability is 
  #                      W(p_m) = -log_M (p_m), m=1, ..., M,
  #     where  W  is a function defined over the same range as ARGVALS.
  #  The fitting criterion is penalized mean squared error:
  #    PENSSE(Wlambda) <- \sum w_i[y_i - f(x_i)]^2 +
  #                     \Wlambda * \int [L W(x)]^2 dx
  #  where L is a linear differential operator defined in argument Lfd,
  #  and w_i is a positive weight applied to the observation.
  #  The function W(x) is expanded by the basis in functional data object
  #    Wfd.
  
  #  Arguments:
  #  ARGVALS ...  Argument value array of length NBIN, the number of 
  #               surprisal values for each curve.  It is assumed that
  #               that these argument values are common to all observed
  #               curves.  
  #  WBIN    ...  A matrix containingg the values to be fit.
  #               This will be an NBIN by M matrix, where NBIN is the 
  #               number of bins containing choice probabilities and M is
  #               the number of options in a specific question or rating
  #               scale.
  #  BMAT0   ...  An initial K by M-1 matrix defining the surprisal curves
  #               via spline functions.  K is the number of basis functions
  #               in the spline expansions, and M is the number of choices
  #               for a particular question in a test or rating scale.
  #  WFDPAR  ...  A functional parameter or fdPar object.  This object
  #               contains the specifications for the functional data
  #               object to be estimated by smoothing the data.  
  #               Note:  WFDPAR is only a container for its 
  #               functional basis object WBASIS, the penalty matrix WPEN, 
  #               and the smoothing parameter Wlambda.  A coefficient
  #               matrix in WFDPAR defined by using a function data object
  #               is discarded, and overwritten by argument BMAT0.
  #  WTVEC   ...  a vector of weights, a vector of N one's by default.
  #  CONV    ...  convergence criterion, 0.0001 by default
  #  ITERLIM ...  maximum number of iterations, 50 by default.
  #  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
  #               no output, if 1, output at each iteration, if higher,
  #               output at each line search iteration. 1 by default.
  #               enabling this option.
  
  #  Returns are:
  #  WFD     ...  Functional data object for W.
  #               Its coefficient matrix an N by NCURVE (by NVAR) matrix
  #               (or array), depending on whether the functional
  #               observations are univariate or multivariate.
  #  FLIST ... A list object or a vector of list objects, one for
  #            each curve (and each variable if functions are multivariate).
  #            Each list object has slots:
  #                 f    ... The sum of squared errors
  #                 grad ... The gradient
  #                 norm ... The norm of the gradient
  #  When multiple curves and variables are analyzed, the lists containing
  #  FLIST objects are indexed linear with curves varying inside
  #  variables.
  
  #  Last modified 18 March 2021 by Jim Ramsay
  
  #  check ARGVALS, a vector of length n
  
  if (!is.numeric(argvals)) stop("ARGVALS is not numeric.")
  argvals <- as.vector(argvals)
  if (length(argvals) < 2) stop("ARGVALS does not contain at least two values.")
  n       <- length(argvals)
  onesobs <- matrix(1,n,1)
  
  #  Check Wbin, an n by M-1 matrix of surprisal values.  
  #  It may not contain negative values.
  
  Wbin <- as.matrix(Wbin)
  Wbindim <- dim(Wbin)
  M <- Wbindim[2]
  if (Wbindim[1] != n) 
      stop("The length of ARGVALS  and the number of rows of WBIN differ.")
  # if (min(Wbin) < 0) stop("WBIN contains negative values.")
  
  #  Check WfdPar and extract WBASIS, WNBASIS, Wlambda and WPENALTY.  
  #  Note that the coefficient matrix is not used.
  
  WfdPar   <- fdParcheck(WfdPar)
  Wbasis   <- WfdPar$fd$basis
  Wnbasis  <- Wbasis$nbasis
  Wlambda  <- WfdPar$lambda
  Wpenalty <- eval.penalty(WfdPar)
  
  #  Check BMAT0, the WNBASIS by M-1 coefficient matrix
  
  if (is.null(Bmat0)) stop("BMAT0 is  NULL.")
  
  Bmatdim <- dim(Bmat0)
  if (Bmatdim[1] != Wnbasis) 
    stop("The first dimension of BMAT0 is not equal to WNBASIS.")
  if (Bmatdim[2] != M-1) 
    stop("The second dimension of BMAT0 is not equal to M - 1.")
  
  #  convert Bmat0 to a vector and NPAR to its length
  
  cvec <- as.vector(Bmat0)
  npar <- length(cvec)
  
  #  Set up the transformation from dimension M-1 to M
  #  where M-vectors sum to zero
  
  M <- dim(Bmat0)[2] + 1
  if (M == 2) {
    root2 <- sqrt(2)
    Zmat <- matrix(1/c(root2,-root2),2,1)
  } else {
    Zmat <- zerobasis(M)
  }
  
  #  Set up the matrix of basis function values 
  
  Phimat <- fda::eval.basis(argvals, Wbasis)
  
  #  check WTVEC
  
  if (is.null(wtvec)) wtvec<-rep(1,n)
  wtvec <- fda::wtcheck(n, wtvec)$wtvec
  
  #  set bounds on coefficient values
  
  climit  <- c(rep(-20,npar),rep(20,npar))
  active  <- 1:Wnbasis
  
  #  initialize matrix Kmat defining penalty term
  
  if (Wlambda > 0) 
  {
    Kmat <- Wlambda*Wpenalty
  } else {
    Kmat <- matrix(0,Wnbasis,Wnbasis)
  }
  
  #  --------------------------------------------------------------------
  #              loop through variables and curves
  #  --------------------------------------------------------------------
  
  #  evaluate log likelihood
  #    and its derivatives with respect to these coefficients
  
  result    <- PENSSEfun(argvals, Wbin, wtvec, Bmat0, Kmat, Zmat, Phimat)
  PENSSE    <- result[[1]]
  DPENSSE   <- result[[2]]
  D2PENSSE  <- result[[3]]
  
  #  compute initial badness of fit measures
  
  f0    <- PENSSE
  gvec0 <- DPENSSE
  hmat0 <- D2PENSSE
  
  # print("gvec0:")
  # print(round(t(gvec0),4))
  # print("hmat0[1:10,1:10]:")
  # print(round(hmat0[1:10,1:10],4))
  
  Foldlist <- list(f = f0, grad = gvec0, norm = sqrt(mean(gvec0^2)))
  
  #  evaluate the initial update vector for correcting the initial bmat
  
  deltac   <- -solve(hmat0,gvec0)
  # print("deltac:")
  # print(round(t(deltac),4))
  cosangle <- -sum(gvec0*deltac)/sqrt(sum(gvec0^2)*sum(deltac^2))
  
  #  initialize iteration status arrays
  
  iternum <- 0
  status <- c(iternum, Foldlist$f, Foldlist$norm)
  if (dbglev >= 1) {
    cat("\n")
    cat("\nIter.   PENSSE   Grad Length")
    cat("\n")
    cat(iternum)
    cat("        ")
    cat(round(status[2],4))
    cat("      ")
    cat(round(status[3],4))
  }
  
  #  -------  Begin iterations  -----------
  
  MAXSTEPITER <- 10
  MAXSTEP     <- 10
  trial       <- 0.2
  linemat     <- matrix(0,3,5)
  gvec        <- gvec0
  dbgwrd      <- dbglev > 1
  
  for (iter in 1:iterlim) {
    iternum <- iternum + 1
    #  take optimal stepsize
    dblwrd <- rep(FALSE,2)
    limwrd <- rep(FALSE,2)
    stpwrd <- FALSE
    ind    <- 0
    Flist  <- Foldlist
    #  compute slope at 0 for line search
    linemat[2,1] <- sum(deltac*Flist$grad)
    #  normalize search direction vector
    sdg     <- sqrt(sum(deltac^2))
    deltac  <- deltac/sdg
    linemat[2,1] <- linemat[2,1]/sdg
    #  break with error condition if initial slope is nonnegative
    if (linemat[2,1] >= 0) {
      if (dbglev > 0) {
        cat("\n")
        print("Initial slope nonnegative.")
        cat("\n")
      }
      ind <- 3
      break
    }
    #  return successfully if initial slope is very small
    if (linemat[2,1] >= -1e-5) {
      if (dbglev > 1) print("Initial slope too small")
      break
    }
    # initialize line search vectors
    linemat[,1:4] <- outer(c(0, linemat[2,1], Flist$f),rep(1,4))
    stepiter <- 0
    if (dbglev >= 2) {
      cat("\n")
      cat(paste("                 ", stepiter, "  "))
      cat(format(round(t(linemat[,1]),4)))
    }
    #  first step set to trial
    ips <- 0
    linemat[1,5]  <- trial
    #  Main iteration loop for linesrch
    for (stepiter in 1:MAXSTEPITER) {
      #  ensure that step does not go beyond limits on parameters
      limflg  <- 0
      #  check the step size
      result <- stepchk(linemat[1,5], cvec, deltac, limwrd, ind,
                        climit, active, dbgwrd)
      linemat[1,5] <- result[[1]]
      ind          <- result[[2]]
      limwrd       <- result[[3]]
      if (linemat[1,5] <= 1e-7) {
        #  Current step size too small  terminate
        Flist  <- Foldlist
        cvecnew <- cvec
        Bmatnew <- matrix(cvecnew,Wnbasis,M-1)
        gvecnew <- gvec
        if (dbglev >= 2) {
          print("Stepsize too small")
          print(linemat[1,5])
        }
        if (limflg) ind <- 1 else ind <- 4
        break
      }
      #  compute new function value and gradient
      cvecnew   <- as.matrix(cvec) + linemat[1,5]*deltac
      Bmatnew   <- matrix(cvecnew,Wnbasis,M-1)
      result    <- PENSSEfun(argvals, Wbin, wtvec, Bmatnew, Kmat, Zmat, Phimat)
      PENSSE    <- result[[1]]
      DPENSSE   <- result[[2]]
      D2PENSSE  <- result[[3]]
      Flist$f   <- PENSSE
      Flist$grad <- DPENSSE
      Flist$norm <- sqrt(mean(DPENSSE^2))
      gvecnew    <- DPENSSE
      # print("gvecnew:")
      # print(round(t(gvecnew),4))
      linemat[3,5] <- Flist$f
      #  compute new directional derivative
      linemat[2,5] <- sum(deltac*gvecnew)
      if (dbglev >= 2) {
        cat("\n")
        cat(paste("                 ", stepiter, "  "))
        cat(format(round(t(linemat[,5]),4)))
      }
      #  compute next step
      result  <- stepit(linemat, ips, dblwrd, MAXSTEP)
      linemat <- result[[1]]
      ips     <- result[[2]]
      ind     <- result[[3]]
      dblwrd  <- result[[4]]
      trial   <- linemat[1,5]
      #  ind == 0 implies convergence
      if (ind == 0 | ind == 5) break
      #  end of line search loop
    }
    cvec <- cvecnew
    Bmat <- matrix(cvec,Wnbasis,M-1)
    gvec <- gvecnew
    status <- c(iternum, Flist$f, Flist$norm)
    if (dbglev > 0) {
      cat("\n")
      cat(iternum)
      cat("        ")
      cat(round(status[2],4))
      cat("      ")
      cat(round(status[3],4))
    }
    #  test for convergence
    if (abs(Flist$f - Foldlist$f) < conv) {
      # cat("\n")
      break
    }
    if (Flist$f >= Foldlist$f) break
    #  compute the Hessian
    hmat <- D2PENSSE
    #  evaluate the update vector
    deltac <- -solve(hmat,gvec)
    cosangle  <- -sum(gvec*deltac)/sqrt(sum(gvec^2)*sum(deltac^2))
    if (cosangle < 0) {
      if (dbglev > 1) print("cos(angle) negative")
      deltac <- -gvec
    }
    Foldlist <- Flist
    #  end of iteration loop
  }
  
  Wfd <- fda::fd(Bmat, Wbasis)
  
  result <- list(Wfd=Wfd, Bmat=Bmat)
  
  return(result)
}

#  ---------------------------------------------------------------

PENSSEfun <- function(argvals, Wbin, wtvec, Bmat, Kmat, Zmat, Phimat) {
  n <- length(argvals)
  M <- dim(Bmat)[2] + 1
  K <- dim(Phimat)[2]
  logM <- log(M)
  onewrd   <- all(wtvec == 1)
  Xmat     <- Phimat %*% Bmat %*% t(Zmat)
  expXmat  <- M^Xmat
  sumexpXmat <- as.matrix(apply(expXmat,1,sum))
  Pmat     <- expXmat/(sumexpXmat %*% matrix(1,1,M))
  Smat     <- -Xmat + (log(sumexpXmat) %*% matrix(1,1,M))/logM
  Rmat     <- Wbin - Smat
  vecBmat  <- matrix(Bmat,K*(M-1),1,byrow=TRUE)
  vecRmat  <- matrix(Rmat,n*M,    1,byrow=TRUE)
  vecKmat  <- kronecker(diag(rep(1,M-1)),Kmat)
  fitscale <- 1
  if (!onewrd) {
    vecwtmat <- diag(rep(wtvec,M))
    PENSSE   <- t(vecRmat) %*% diag(as.numeric(wtvec)) %*% vecRmat/fitscale + 
      t(vecBmat) %*% vecKmat %*% vecBmat
  } else {
    PENSSE   <- t(vecRmat) %*% vecRmat/fitscale + t(vecBmat) %*% vecKmat %*% vecBmat
  }
  DvecXmatDvecB <- kronecker(Zmat,Phimat)
  DvecSmatDvecX <- matrix(0,n*M,n*M)
  m2 <- 0
  for (m in 1:M) {
    m1 <- m2 + 1
    m2 <- m2 + n
    m4 <- 0
    for (l in 1:M) {
      m3 <- m4 + 1
      m4 <- m4 + n
      diagPl <- diag(Pmat[,l])
      DvecSmatDvecX[m1:m2,m3:m4] <- diagPl
    }
  }
  DvecSmatDvecX <- DvecSmatDvecX - diag(rep(1,n*M))
  DvecSmatDvecB <- DvecSmatDvecX %*% DvecXmatDvecB
  if (!onewrd) {
    DPENSSE  <- -2*t(DvecSmatDvecB) %*% vecwtmat %*% vecRmat/fitscale + 
      2*vecKmat %*% vecBmat
  } else {
    DPENSSE  <- -2*t(DvecSmatDvecB) %*% vecRmat/fitscale + 2*vecKmat %*% vecBmat
  }
  if (!onewrd) {
    D2PENSSE <-  
      2*((DvecSmatDvecB) %*% vecwtmat %*% DvecSmatDvecB)/fitscale + 2*vecKmat
  } else {
    D2PENSSE <-  2*(t(DvecSmatDvecB) %*% DvecSmatDvecB)/fitscale  + 2*vecKmat
  }
  
  return(list(
    PENSSE   = PENSSE, 
    DPENSSE  = DPENSSE, 
    D2PENSSE = D2PENSSE)
  )
  
}

# ------------------------------------------------------------------
ycheck <- function(y, n) {
  
  #  check Y
  
  if (is.vector(y)) y <- as.matrix(y)
  
  if (!inherits(y, "matrix") && !inherits(y, "array"))
    stop("Y is not of class matrix or class array.")
  
  ydim <- dim(y)
  
  if (ydim[1] != n) stop("Y is not the same length as ARGVALS.")
  
  #  set number of curves and number of variables
  
  ndim  <- length(ydim)
  if (ndim == 2) {
    ncurve <- ydim[2]
    nvar   <- 1
  }
  if (ndim == 3) {
    ncurve <- ydim[2]
    nvar   <- ydim[3]
  }
  if (ndim > 3) stop("Second argument must not have more than 3 dimensions")
  
  
  return(list(y=y, ncurve=ncurve, nvar=nvar, ndim=ndim))
  
}

# ------------------------------------------------------------------
fdParcheck <- function (fdPar) {
  if (!inherits(fdPar, "fdPar")) {
    if (inherits(fdPar, "fd") || inherits(fdPar, "basisfd")) {
      fdPar <- fdPar(fdPar)
    } else
      stop(paste("'fdPar' is not a functional parameter object,",
                 "not a functional data object, and",
                 "not a basis object."))
  }
  
  return(fdPar)
  
}

# ------------------------------------------------------------------
stepchk <- function(oldstep, cvec, deltac, limwrd, ind,
                    climit=50*c(-rep(1,ncvec), rep(1,ncvec)),
                    active=1:ncvec, dbgwrd) {
  #  stepcheck checks the step size to keep parameters within boundaries
  
  # Last changed 2018 by Jim Ramsay
  
  # define vectors containing lower and upper limits
  
  ncvec   <- length(deltac)
  bot     <- climit[1:ncvec]
  top     <- climit[ncvec+(1:ncvec)]
  
  newstep <- oldstep
  
  #  ensure that step does not go beyond lower limit on parameters
  #  limwrd[2] flags that the lower limit has been hit once
  
  stepi   <- oldstep*deltac
  stepmin <- min(stepi)
  index   <- stepi[active] == stepmin
  if (any(stepi[index] < bot[index]-cvec[index]) &
      any(deltac[index] != 0) )  {
    stepnew <- min((bot[index]-cvec[index])/deltac[index])
    if (dbgwrd) {
      print("Lower limit reached ... new step:")
      cat(c(stepi, round(c(oldstep, stepnew),4)),"\n")
      cat(round(cvec + stepnew*deltac,4),"\n")
    }
    newstep <- stepnew
    if (limwrd[2]) ind <- 1 else limwrd[2] <- TRUE
  } else {
    limwrd[2] <- FALSE
  }
  
  #  check whether upper limit has been reached twice in a row
  
  #  ensure that step does not go beyond upper limit on parameters
  #  limwrd[1] flags that the upper limit has been hit once
  
  stepi   <- oldstep*deltac
  stepmax <- max(stepi)
  index   <- stepi[active] == stepmax
  if (any(stepi[index] > top[index]-cvec[index]) &
      any(deltac[index] != 0) ) {
    stepnew <- min((top[index]-cvec[index])/deltac[index])
    if (dbgwrd) {
      print("Upper limit reached ... new step:")
      cat(c(stepi, round(c(oldstep, stepnew),4)),"\n")
    }
    newstep <- stepnew
    if (limwrd[1]) ind <- 1 else limwrd[1] <- TRUE
  } else {
    limwrd[1] <- FALSE
  }
  
  return(list(newstep, ind, limwrd))
}

# ------------------------------------------------------------------
stepit <- function(linemat, ips, dblwrd, MAXSTEP) {
  #STEPIT computes next step size in line search algorithm
  
  #  Arguments:
  #  LINEMAT:  Row 1 contains step values
  #            Row 2 contains slope values
  #            Row 3 contains function values
  #  IPS:      If 1, previous slope was positive
  #  DBLWRD:   Vector of length 2:  dblwrd[1] TRUE means step halved
  #                                 dblwrd[2] TRUE means step doubled
  #  MAXSTEP:  maximum size of step
  
  #  Last modified 29 June 2018 by Jim Ramsay
  
  #  Wolfe condition 1
  test1.1 <- linemat[3,5] <= linemat[3,1] + linemat[1,5]*linemat[2,1]/20
  #  Wolfe condition 2
  test1.2 <- abs(linemat[2,5]) <= abs(linemat[2,1])/10 
  # test1 <- test1.1 && test1.2
  test1 <- test1.2
  test2 <- linemat[3,5] > linemat[3,1]
  test3 <- linemat[2,5] > 0
  if ((test1 || !test3) && test2) {
    #  ************************************************************
    #  function is worse and either slope is satisfory or negative
    ips <- 0        #  step is halved
    if (dblwrd[2]) {
      ind <- 5
      return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
    }
    linemat[1,5] <- min(c(linemat[1,5]/2, MAXSTEP))
    linemat[,2] <- linemat[,1]
    linemat[,3] <- linemat[,1]
    dblwrd[1] <- TRUE
    ind <- 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  *********************************************************
  if (test1) {
    #  test1 means successful convergence
    ind <- 0
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  **********************************************************
  if (test3) {
    #  Current slope is positive
    ips <- 1
    linemat[,4] <- linemat[,5]
    deltaf <- linemat[3,3] - linemat[3,5]
    z <- (3/(linemat[1,5] - linemat[1,3]))*deltaf + linemat[2,3] + linemat[2,5]
    w <- z * z - linemat[2,3] * linemat[2,5]
    if (abs(linemat[2,3] + linemat[2,5] + 2 * z) >= 1e-05 && w > 0) {
      w <- sqrt(w)
      linemat[1,5] <- linemat[1,3] + (1 - ((linemat[2,5] + w - z)/ 
                                            (linemat[2,5] - linemat[2,3] + 2 * w))) * (linemat[1,5] - linemat[1,3])
    } else {
      #  linear interpolation necessary
      aerror <- linemat[1,3]
      if (linemat[1,5] > linemat[1,3]) {
        aerror <- linemat[1,5]
      }
      linemat[1,5] <- linemat[1,3] - linemat[2,3] * 
        ((linemat[1,5] - linemat[1,3])/ 
           (linemat[2,5] - linemat[2,3]))
      if (linemat[1,5] > 2 * aerror) {
        linemat[1,5] <- 2 * aerror
      }
    }
    linemat[1,5] <- min(c(linemat[1,5], MAXSTEP))
    dblwrd <- rep(FALSE,2)
    ind <- 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  *************************************************************
  #  Current slope is negative or zero
  linemat[,2] <- linemat[,3]
  linemat[,3] <- linemat[,5]
  if (ips == 1) {
    #  *****************************************************
    #  previous slope is positive
    deltaf <- linemat[3,5] - linemat[3,4]
    z <- c(3/(linemat[1,4] - linemat[1,5])) * deltaf + 
      linemat[2,5] + linemat[2,4]
    w <- z * z - linemat[2,5] * linemat[2,4]
    if (abs(linemat[2,5] + linemat[2,4] + 2 * z) >= 1e-05 && w > 0) {
      w <- sqrt(w)
      linemat[1,5] <- linemat[1,5] + (1 - ((linemat[2,4] + w - z)/ 
                                            (linemat[2,4] - linemat[2,5] + 
                                               2 * w))) * (linemat[1,4] - linemat[1,5])
    } else {
      aerror <- linemat[1,5]
      if (linemat[1,4] > linemat[1,5]) {
        aerror <- linemat[1,4]
      }
      linemat[1,5] <- linemat[1,5] - linemat[2,5] * 
        ((linemat[1,4] - linemat[1,5])/ 
           (linemat[2,4] - linemat[2,5]))
      if (linemat[1,5] > 2 * aerror) {
        linemat[1,5] <- 2 * aerror
      }
    }
    linemat[1,5] <- min(c(linemat[1,5], MAXSTEP))
    dblwrd <- rep(FALSE,2)
    ind <- 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  ******************************************************
  if ((linemat[2,3] - linemat[2,2]) * (linemat[1,3] - linemat[1,2]) > 0) {
    #  previous slope is negative
    z <- c(3/(linemat[1,3] - linemat[1,2])) * (linemat[3,2] - linemat[3,3]) + 
      linemat[2,2] + linemat[2,3]
    w <- z * z - linemat[2,2] * linemat[2,3]
    if (abs(linemat[2,2] + linemat[2,3] + 2 * z) >= 1e-05 && w > 0) {
      w <- sqrt(w)
      linemat[1,5] <- linemat[1,2] + 
        (1 - ((linemat[2,3] + w - z)/(linemat[2,3] - linemat[2,2] + 
                                        2 * w))) * (linemat[1,3] - linemat[1,2])
    } else {
      linemat[1,5] <- linemat[1,2] - linemat[2,2] * 
        ((linemat[1,3] - linemat[1,2])/ 
           (linemat[2,3] - linemat[2,2]))
    }
    linemat[1,5] <- min(c(linemat[1,5], MAXSTEP))
    dblwrd <- rep(FALSE,2)
    ind <- 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  } else {
    #  previous slope also negative but not as much
    if (dblwrd[1]) {
      ind <- 5
      return(
        list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
    } else {
      linemat[1,5] <- 2 * linemat[1,5]
      linemat[1,5] <- min(c(linemat[1,5], MAXSTEP))
      dblwrd[2] <- TRUE
      ind <- 2
      return(
        list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
    }
  }
  ind <- 2
  
}

# ------------------------------------------------------------------

fdParcheck = function (fdParobj) {
  if (!inherits(fdParobj, "fdPar")) {
    if (inherits(fdParobj, "fd") || inherits(fdParobj, "basisfd")) {
      fdParobj <- fdPar(fdParobj)
    } else
      stop(paste("'fdParobj' is not a functional parameter object,",
                 "not a functional data object, and",
                 "not a basis object."))
  }
  return(fdParobj)
}

# ------------------------------------------------------------------

zerobasis <- function(k) {
# ZEROBASIS constructes a K by K-1 matrix that maps an unrestricted matrix B with K - 1 rows by 
#  the linear transformation ZEROBASIS %*% B = C into the subspace of matrices with K rows having #  column sums equal to zero.  
#  The matrix has orthonormal columns, so that crossprod(ZEROBASIS) is the identity matrix
#  of order K - 1.

  tk <- 0:(k-1) + 0.5
  fbasis     <- fda::create.fourier.basis(k,k)
  fbasmat    <- fda::eval.basis(tk, fbasis)
  fbasmat    <- fbasmat[,2:k]
  fbasnorm   <- sqrt(apply(fbasmat^2,2,sum))
  zerobasmat <- fbasmat/outer(rep(1,k),fbasnorm)
  return(zerobasmat)
}
