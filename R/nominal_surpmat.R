nominal_surpmat <- function(parmat, bdry=c(-4,4), nderiv=0) {
  # Inputs mirt parameter matrix and computes values of probability
  # and surprisal over a fine mesh of values of theta
  # Last modified 21 April 2023
  
  M       <- nrow(parmat)
  logM    <- log(M)
  nfine   <- 101
  indfine <- seq(bdry[1],bdry[2],len=nfine)
  expmat  <- matrix(0,nfine,M)
  for (m in 1:M) {
    linmat     <- parmat[m,1]*indfine + parmat[m,2]
    expmat[,m] <- exp(linmat)
  }
  prbmat  <- expmat/apply(expmat,1,sum)
  surpmat <- -log(prbmat)/logM
  if (nderiv == 0) {
    return(list(prbmat=prbmat, surpmat=surpmat, expmat=expmat,
                infovec=NULL, Dbsurparray=NULL, Dbsurparray=NULL,
                Dthetasurpmat=NULL, D2thetasurpmat=NULL))
  } else {
    eyeM <- diag(rep(1,M))
    oneM <- rep(1,M)
    Dasurparray <- array(0, c(nfine,M,M))
    Dbsurparray <- array(0, c(nfine,M,M))
    Dsurpmatj   <- matrix(0,M,M)
    slpmat      <- outer(rep(1,nfine),parmat[,1])
    denvec      <- apply(expmat,1,sum)
    denmat      <- outer(denvec,oneM)
    for (j in 1:nfine) {
      Dsurpmatj        <- 
        (-eyeM + diag(expmat[j,])/outer(oneM,denmat[j,]))/logM
      Dasurparray[j,,] <- indfine[j]*Dsurpmatj
      Dbsurparray[j,,] <-          Dsurpmatj
    }
    Dthetasurpmat  <- (-slpmat + (slpmat*expmat)/denmat)/logM
    D2thetasurpmat <- 
      ((slpmat*expmat)/denmat - (slpmat*expmat*denmat)/denmat^2)/logM
    infovec        <- pracma::cumtrapz(indfine, sqrt(apply(Dthetasurpmat^2,1,sum)))
    return(list(infovec=infovec, prbmat=prbmat, surpmat=surpmat, expmat=expmat,
                Dasurparray=Dasurparray, Dbsurparray=Dbsurparray,
                Dthetasurpmat=Dthetasurpmat, D2thetasurpmat=D2thetasurpmat))
  }
}
