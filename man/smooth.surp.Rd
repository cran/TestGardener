\name{smooth.surp}
\alias{smooth.surp}
\title{
  Fit data with surprisal smoothing.
}
\description{
  Surprisal is -log(probability) where the logarithm is to the base being the 
  dimension \code{M} of the multinomial observation vector.  The surprisal 
  curves for each question are estimated by fitting the surprisal values of 
  binned data using curves whose values are within the \code{M-1}-dimensional 
  surprisal subspace that is within the space of non-negative M-dimensional 
  vectors.
}
\usage{
smooth.surp(binctr, Sbin, Bmat, Sbasis, Zmat, wtvec=NULL, conv=1e-4,
            iterlim=50, dbglev=0)
}
\arguments{
  \item{binctr}{
    Argument value array of length N, where N is the number of observed
    curve values for each curve.  It is assumed that that these argument
    values are common to all observed curves.  If this is not the case,
    you will need to run this function inside one or more loops,
    smoothing each curve separately.
  }
  \item{Sbin}{
    A \code{nbin} by \code{M_i} matrix of surprisal values to be fit.  
  }
  \item{Bmat}{A \code{Snbasis} by \code{M_i - 1} matrix containing starting
  values for the iterative optimization of the least squares fit of the 
  surprisal curves to the surprisal data.
  }
  \item{Sbasis}{
    A functional data basis object.  
  }
  \item{Zmat}{An \code{M by M-1} matrix satisfying \code{Zmat'Zmat <- I} and 
    \code{Zmat'1 <- 0}.}
  \item{wtvec}{
    A vector of weights to be used in the smoothing.
  }
  \item{conv}{
    A convergence criterion.
  }
  \item{iterlim}{
    the maximum number of iterations allowed in the minimization of
    error sum of squares.
  }
  \item{dbglev}{
    Either 0, 1, or 2.  This controls the amount information printed out
    on each iteration, with 0 implying no output, 1 intermediate output
    level, and 2 full output.  If either level 1 or 2 is specified, it
    can be helpful to turn off the output buffering feature of S-PLUS.
  }
}
\value{
  A named list of class \code{surpFd} with these members:
  \item{PENSSE}{The final value of the penalized fitting criterion.}
  \item{DPENSSE}{The final gradient of the penalized fitting criterion.}
  \item{D2PENSSE}{The final hessian of the fitting criterion.}
  \item{SSE}{The final value of the error sum of squares.}
  \item{DSSE}{The final gradient of the error sum of squares.}
  \item{D2SSE}{The final hessian of the error sum of squares.}
  \item{DvecSmatDvecB}{The final cross derivative DvecSmatDvecX times 
    DvecXmatDvecB of the surprisal curve and the basis coordinates.}
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{eval.surp}},
  \code{\link{ICC_plot}},
  \code{\link{Sbinsmth}}
}
\examples{
  oldpar <- par(no.readonly=TRUE)
  # Assemble the objects in the list arguments Bmat and surpList
  SfdList1 <- Quant_13B_problem_parmList$SfdList[[1]]
  binctr   <- Quant_13B_problem_parmList$binctr
  M        <- SfdList1$M
  Bmat     <- SfdList1$Sfd$coef
  Sbasis   <- SfdList1$Sfd$basis
  Snbasis  <- Sbasis$nbasis
  Phimat   <- fda::eval.basis(binctr, Sbasis)
  Zmat     <- SfdList1$Zmat
  Sbin     <- SfdList1$Sbin
  # add some noise by rounding
  Bmat     <- round(Bmat,0)
  Kmat     <- matrix(0,Snbasis,Snbasis)
  wtvec    <- NULL
  surpList <- list(binctr=binctr, Sbin=Sbin, wtvec=wtvec, 
                 Kmat=Kmat, Zmat=Zmat, Phimat=Phimat, M=M)
  Bvec     <- matrix(Bmat,Sbasis$nbasis*(M-1),1,byrow=TRUE)
  #  run surp.fit to get initial values
  result   <- surp.fit(Bvec, surpList)
  print(paste("Initial error sum of squares =",round(result$SSE,3)))
  print(paste("Initial gradient norm =",round(norm(result$DSSE),5)))
  #  optimize SSE
  result   <- smooth.surp(binctr, Sbin, Bmat, Sbasis, Zmat)
  print(paste("Optimal error sum of squares =",round(result$SSE,3)))
  print(paste("Optimal gradient norm =",round(norm(result$DSSE),5)))
  par(oldpar)
}

