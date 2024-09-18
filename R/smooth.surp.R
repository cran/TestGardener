smooth.surp <- function(binctr, Sbin, Bmat, Sbasis, Zmat, wtvec=NULL, 
                        conv=1e-4, iterlim=50, dbglev=0) {
  #  Smooths the relationship of SBIN to BINCTR within a single item
  #     using weights in STVEC by fitting surprisal functions 
  #     to a set of surprisal transforms of choice probabilities, 
  #     where the surprisal transformation of each probability is 
  #             S(p_m) = -log_M (p_m), m=1, ..., M,
  #     where  S  is a function defined over the same range as binctr.
  #  The fitting criterion is penalized mean squared error:
  #              PENSSE <- \sum w_i[y_i - f(x_i)]^2
  #  where L is a linear differential operator defined in argument Lfd,
  #  and w_i is a positive weight applied to the observation.
  #
  #  This version uses Numerical Recipes function lnsrch()
  #
  #  Arguments:
  #  BINCTR  ...  Argument value array of length NBIN, the number of 
  #               surprisal values for each curve.  It is assumed that
  #               that these argument values are common to all observed
  #               curves.  
  #  SBIN    ...  A matrix containing the values to be fit.
  #               This will be an NBIN by M matrix, where NBIN is the 
  #               number of bins containing choice probabilities and M is
  #               the number of options in a specific question or rating
  #               scale.
  #  BMAT   ...  An initial K by M-1 matrix defining the surprisal curves
  #               via spline functions.  K is the number of basis functions
  #               in the spline expansions, and M is the number of choices
  #               for a particular question in a test or rating scale.
  #  SBASIS ...  B-spline basis object for representing the surprisal curves.
  #  ZMAT    ... An M by M-1 matrix satisfying Z'Z = I and Z'1 = 0.
  #  WTVEC   ...  a vector of weights, a vector of N one's by default.
  #  CONV    ...  convergence criterion, 0.0001 by default
  #  ITERLIM ...  maximum number of iterations, 50 by default.
  #  DBGLEV  ...  Controls the level of output on each iteration.  If 0,
  #               no output, if 1, output at each iteration, if higher,
  #               output at each line search iteration. 1 by default.
  #               enabling this option.
  
  #  Returns are:
  #  FLIST ... A list object or a vector of list objects, one for
  #            each curve (and each variable if functions are multivariate).
  #            Each list object has slots:
  #                 f    ... The sum of squared errors
  #                 grad ... The gradient
  #                 norm ... The norm of the gradient
  #  When multiple curves and variables are analyzed, the lists containing
  #  FLIST objects are indexed linear with curves varying inside
  #  variables.
  
  #  Last modified 29 March 2024 by Jim Ramsay
  
  
  
  #-----------------------------------------------------------------------------
  #                             check arguments   
  #-----------------------------------------------------------------------------
  
  if (!is.numeric(binctr)) stop("binctr is not numeric.")
  binctr <- as.vector(binctr)
  if (length(binctr) < 2) stop("binctr does not contain at least two values.")
  nbin    <- length(binctr)
  onesobs <- matrix(1,nbin,1)
  
  #  Check Sbin, an nbin by M-1 matrix of surprisal values.  
  #  It may not contain negative values.
  
  Sbin    <- as.matrix(Sbin)
  Sbindim <- dim(Sbin)
  M       <- Sbindim[2]
  Pbin    <- M^(-Sbin)
  if (Sbindim[1] != nbin) 
      stop("The length of binctr  and the number of rows of Sbin differ.")
  Snbasis  <- Sbasis$nbasis
  
  #  Check Bmat, the SNBASIS by M-1 coefficient matrix
  
  if (is.null(Bmat)) stop("Bmat is  NULL.")
  
  Bmatdim <- dim(Bmat)
  if (Bmatdim[1] != Snbasis) 
    stop("The first dimension of Bmat is not equal to SNBASIS.")
  if (Bmatdim[2] != M-1) 
    stop("The second dimension of Bmat is not equal to M - 1.")
  
  #  convert Bmat to a vector and NPAR to its length
  
  K    <- Snbasis
  Bvec <- matrix(Bmat,K*(M-1),1,byrow=TRUE)
  npar <- length(Bvec)
  
  #  Set up the matrix of basis function values 
  
  Phimat <- fda::eval.basis(binctr, Sbasis)
  
  #  check WTVEC
  
  if (is.null(wtvec)) wtvec<-rep(1,nbin)
  wtvec <- fda::wtcheck(nbin, wtvec)$wtvec
  
  #  initialize matrix Kmat defining penalty term
  
  Kmat <- matrix(0,Snbasis,Snbasis)
  
  #-----------------------------------------------------------------------------
  #  Set up initial list object for data required by function surp.fit()
  #-----------------------------------------------------------------------------
  
  surpList <- list(binctr=binctr, Sbin=Sbin, wtvec=wtvec, 
                   Kmat=Kmat, Zmat=Zmat, Phimat=Phimat, M=M)
  
  #  evaluate log likelihood and its derivatives with respect to these 
  #  coefficients and compute initial badness of fit measures
  
  Bvecold <- matrix(Bmat, Snbasis*(M-1),1)
  #  ----------------------------------
  result <- surp.fit(Bvecold, surpList)
  #  ----------------------------------
  
  fold    <- result$PENSSE
  gvec    <- result$DPENSSE
  hmat    <- result$D2PENSSE
  
  Flist    <- list(f = fold, grad = gvec, norm = sqrt(mean(gvec^2)))
  Foldlist <- Flist
  
  #  evaluate the initial update vector for correcting the initial bmat
  
  pvec     <- -solve(hmat,gvec)
  cosangle <- -sum(gvec*pvec)/sqrt(sum(gvec^2)*sum(pvec^2))
  
  #  initialize iteration status arrays
  
  iternum <- 0
  status <- c(iternum, Foldlist$f, Foldlist$norm)
  if (dbglev > 0) {
    cat("\n")
    cat("\nIter.   PENSSE   Grad Length")
    cat("\n")
    cat(iternum)
    cat("        ")
    cat(round(status[2],4))
    cat("      ")
    cat(round(status[3],4))
  }
  
  #-----------------------------------------------------------------------------
  #                            Begin iterations  
  #-----------------------------------------------------------------------------
    
  STEPMAX <- 1
  iternum <- 0
  for (iter in 1:iterlim) {
    iternum <- iternum + 1
    if (any(is.na(pvec))) {
      FList   <- Foldlist
      Bvec    <- Bvecold
      #  ----------------------------------
      result <- surp.fit(Bvec, surpList)
      #  ----------------------------------
      f        <- result$PENSSE
      gvec     <- result$DPENSSE 
      hmat     <- result$D2PENSSE
      SSE      <- result$SSE
      DSSE     <- result$DSSE
      D2SSE    <- result$D2SSE
      Flist$f    <- f
      Flist$grad <- gvec
      Flist$norm <- sqrt(mean(gvec^2))
      break
    } else {
      #  take optimal stepsize
      gvecnorm <- sqrt(sum(gvec^2))
      # print(gvecnorm)
      if (gvecnorm > conv) {
        lnsrch_result <- 
          fda::lnsrch(Bvecold, fold, gvec, pvec, surp.fit, surpList, STEPMAX)
        check    <- lnsrch_result$check
        if (check) stop("lnsrch failure")
      } else {
        Bvec     <- lnsrch_result$x
        gvecnorm <- norm(as.numeric(lnsrch_result$gvec))
        print(gvecnorm)
        print("lnsrch failure")
        print("Bvecold")
        print(as.numeric(lnsrch_result$Bvec))
        print(paste("fold =", as.numeric(lnsrch_result$f)))
        print("gvec")
        print(as.numeric(lnsrch_result$gvec))
        print("pvec")
        print(as.numeric(lnsrch_result$pvec))
        print(paste("slope(gvec.pvec =", as.numeric(sum(gvec*pvec))))
        break
      }
      #  set up new Bvec and evaluate function, gradient and hessian
      #  ----------------------------------
      result <- surp.fit(Bvec, surpList)
      #  ----------------------------------
      f        <- result$PENSSE
      gvec     <- result$DPENSSE 
      hmat     <- result$D2PENSSE
      SSE      <- result$SSE
      DSSE     <- result$DSSE
      D2SSE    <- result$D2SSE
      #  set up list object for current fit
      Flist$f    <- f
      Flist$grad <- gvec
      Flist$norm <- sqrt(mean(gvec^2))
      #  store current values for next iteration
      Bvecold <- Bvec
      fold    <- f
      #  display results at this iteration if dbglev > 0
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
        break
      }
      #  also terminate iterations if new fit is worse than old
      if (Flist$f >= Foldlist$f) break
      #  set up objects for new iteration
      #  evaluate the update vector using Newton Raphson direction
      pvec      <- -solve(hmat,gvec)
      cosangle  <- -sum(gvec*pvec)/sqrt(sum(gvec^2)*sum(pvec^2))
      if (cosangle < 0) {
        if (dbglev > 1) print("cos(angle) negative")
        pvec <- -gvec
      }
      Foldlist <- Flist
    }
  }
  if (dbglev > 0) cat("\n")

  #  end of iteration loop, output results
  
  Bmat <- matrix(Bvec, Snbasis, M-1)
  # Sfd  <- fda::fd(Bmat, Sbasis)
  #  -------------------------------
  result <- surp.fit(Bvec, surpList)
  #  -------------------------------
  f             <- result$PENSSE
  gvec          <- result$DPENSSE 
  hmat          <- result$D2PENSSE
  PENSSE        <- result$PENSSE
  DPENSSE       <- result$DPENSSE 
  D2PENSSE      <- result$D2PENSSE
  SSE           <- result$SSE
  DSSE          <- result$DSSE
  D2SSE         <- result$D2SSE
  DvecSmatDvecB <- result$DvecSmatDvecB
  RMSE          <- result$RMSE
  Entropy       <- result$Entropy
  
  surpFd <- list(Bmat=Bmat, f=f, gvec=gvec, hmat=hmat,
                 PENSSE=PENSSE, DPENSSE=DPENSSE, D2PENSSE=D2PENSSE,
                 SSE=SSE, DSSE=DSSE, D2SSE=D2SSE,
                 DvecSmatDvecB=DvecSmatDvecB, RMSE=RMSE, Entropy=Entropy)
  class(surpFd) <- 'surpfd'
  
  return(surpFd)
}


