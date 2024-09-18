surp.fit <- function(Bvec, surpList) {
  
  #  This function is called within smooth.surp() to
  #  evaluate fit within a single item at each iteration
  
  # Arguments:
  #. Bvec     ... A K by M-1 matrix Bmat of coefficients for surprisal curves
  #               in row-wise column vector format for use with 
  #               with function lnsrch().
  #. surpList ... a list containing objects required to assess fit: 
  #               M, binctr, Sbin, wtvec, Kmat, Zmat and  Phimat.
  
  #. Last modified 27 March 2024 by Jim

  #  extract objects from surpList
  
  M       <- surpList$M
  binctr  <- surpList$binctr 
  Sbin    <- surpList$Sbin 
  wtvec   <- surpList$wtvec 
  Kmat    <- surpList$Kmat
  Zmat    <- surpList$Zmat 
  Phimat  <- surpList$Phimat 
  
  # set up dimensions
  
  nbin    <- length(binctr)
  K       <- dim(Phimat)[2]
  Pbin    <- M^(-Sbin)
  vecPbin <- matrix(Pbin,nbin*M,1,byrow=TRUE)
  
  #  compute fit, gradient and hessian
  
  logM     <- log(M)
  onewrd   <- all(wtvec == 1)
  Bmat     <- matrix(Bvec,K,M-1)
  Xmat     <- Phimat %*% Bmat %*% t(Zmat)
  expXmat  <- M^Xmat
  sumexpXmat <- as.matrix(apply(expXmat,1,sum))
  Pmat     <- expXmat/(sumexpXmat %*% matrix(1,1,M))
  Smat     <- -Xmat + (log(sumexpXmat) %*% matrix(1,1,M))/logM
  Rmat     <- Sbin - Smat
  RMSE     <- apply(Pmat*(Rmat^2),1,mean)
  Entropy  <- apply(Pmat*Rmat,1,sum)
  vecRmat  <- matrix(Rmat,nbin*M,    1,byrow=TRUE)
  vecKmat  <- kronecker(diag(rep(1,M-1)),Kmat)
  fitscale <- 1
  #  compute raw fit and its penalized version
  vecwtmat <- diag(matrix(Pbin,nbin*M,1,byrow=TRUE))
  SSE      <- crossprod(vecRmat*vecPbin,vecRmat)
  PENSSE   <- SSE/fitscale + t(Bvec) %*% vecKmat %*% Bvec
  DvecXmatDvecB <- kronecker(Zmat,Phimat)
  DvecSmatDvecX <- matrix(0,nbin*M,nbin*M)
  m2 <- 0
  for (m in 1:M) {
    m1 <- m2 + 1
    m2 <- m2 + nbin
    m4 <- 0
    for (l in 1:M) {
      m3 <- m4 + 1
      m4 <- m4 + nbin
      diagPl <- diag(Pmat[,l])
      DvecSmatDvecX[m1:m2,m3:m4] <- diagPl
    }
  }
  DvecSmatDvecX <- DvecSmatDvecX - diag(rep(1,nbin*M))
  DvecSmatDvecB <- DvecSmatDvecX %*% DvecXmatDvecB
  DSSE     <- -2*t(DvecSmatDvecB) %*% vecRmat
  DPENSSE  <- DSSE/fitscale + 2*vecKmat %*% Bvec
  D2SSE    <- 2*crossprod(DvecSmatDvecB)
  D2PENSSE <- D2SSE/fitscale + 2*vecKmat
  
  #  return list object containing raw and penalized fit data
  
  return(
    list(
      PENSSE = PENSSE, DPENSSE = DPENSSE, D2PENSSE = D2PENSSE,
      SSE    = SSE,    DSSE    = DSSE,    D2SSE    = D2SSE,
      Rmat   = Rmat,   RMSE    = RMSE,    Entropy  = Entropy,
      DvecSmatDvecB = DvecSmatDvecB)
    )
}

