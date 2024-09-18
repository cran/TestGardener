\name{surp.fit}
\alias{surp.fit}
\title{
Objects resulting for assessing fit of surprisal matrix to surprisal data
}
\description{
  This function is called by function \code{smooth.surp()} and computes the
  penalized verstion of the objective function value, its derivative vector 
  and the second derivative matrix, as well as their unpenalized versions.
  Also returned are alternative fitting objects: the residual matrix,
  the root-mean-square of the matrix fit, and the entropy value.
}
\usage{
surp.fit(Bvec, surpList)
}
\arguments{
  \item{Bvec}{The K by M-1 parameter matrix defining the fit to the data
    in row-wise column vector format for use with function lnsrch().}
  \item{surpList}{A list object containing objects M, binctr, Sbin, wtvec, 
    Kmat, Zmat and  Phimat.}
}
\value{
  A named list of class \code{surpFd} with these members:
  \item{PENSSE}{value of the penalized fitting criterion.}
  \item{DPENSSE}{gradient of the penalized fitting criterion.}
  \item{D2PENSSE}{hessian of the fitting criterion.}
  \item{SSE}{value of the error sum of squares.}
  \item{DSSE}{gradient of the error sum of squares.}
  \item{D2SSE}{hessian of the error sum of squares.}
  \item{DvecSmatDvecB}{cross derivative DvecSmatDvecX times 
    DvecXmatDvecB of the surprisal curve and the basis coordinates.}
  \item{Rmat}{residual matrix for the fit to the surprisal matrix.}
  \item{RMSE}{root-mean-squared scalar fit value.}
  \item{Entropy}{entropy of the fit to the data.}
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
  \code{\link{smooth.surp}},
  \code{\link{Sbinsmth}}
}
\examples{
# Assemble the objects in the list arguments Bmat and surpList
SfdList1 <- Quant_13B_problem_parmList$SfdList[[1]]
Bmat     <- SfdList1$Sfd$coef
binctr   <- Quant_13B_problem_parmList$binctr
M        <- SfdList1$M
Sbasis   <- SfdList1$Sfd$basis
Zmat     <- SfdList1$Zmat
Sbin     <- SfdList1$Sbin
Phimat   <- fda::eval.basis(binctr, Sbasis)
Snbasis  <- Sbasis$nbasis
Kmat     <- matrix(0,Snbasis,Snbasis)
wtvec    <- NULL
Bvec     <- matrix(Bmat, Snbasis*(M-1),1,byrow=TRUE)
#  display coefficient matrix
print(round(Bmat,2))
#  set up argument surpList
surpList <- list(binctr=binctr, Sbin=Sbin, wtvec=wtvec, 
                 Kmat=Kmat, Zmat=Zmat, Phimat=Phimat, M=M)
#  run surp.fit
result <- surp.fit(Bvec, surpList)
print(paste("Error sum of squares =",round(result$SSE,3)))
print(paste("Gradient norm =",round(norm(result$DSSE),5)))
print("Entropy of item = at bin centres")
print(round(result$Entropy,3))
}