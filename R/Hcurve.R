Hcurve <- function(WfdList, U, nderiv=0) {
  #  HFUN computes the negative log likelihoods for (a set of examinees
  #  at a single value theta.
  #  Items can be either binary or multi-option
  #  The analysis is within [0,100], and does not transform to the real line
  
  #  Last updated 8 November 2022
  
  N <- nrow(U)
  n <- ncol(U)
  
  #  loop through items to compute negative log likelihood values in H
  
  Hfine <- matrix(0,101,N)
  for (j in 1:N) {
    for (item in 1:n) {
      Uveci <- U[,item]
      if (!is.null(Uveci)) {
        WListi <- WfdList[[item]]
        if (nderiv == 0) {
          Wmatfinei <- WListi$Wmatfine
          Hfine[,j] <- Hfine[,j] +   Wmatfinei[,U[j,item]]
        }
        if (nderiv == 1) {
          DWmatfinei <- WListi$DWmatfine
          Hfine[,j] <- Hfine[,j] +  DWmatfinei[,U[j,item]]
        }
        if (nderiv == 2) {
          D2Wmatfinei <- WListi$D2Wmatfine
          Hfine[,j] <- Hfine[,j] + D2Wmatfinei[,U[j,item]]
        }
      }
    }
  }
  
  return(Hfine)
  
}
