\name{entropies}
\alias{entropies}
\title{
Entropy measures of inter-item dependency
}
\description{
Entropy $I_1$ is a scalar measure of how much informationn is required to predict 
the outcome of a choice number 1 exactly, and consequently is a measure of item effectiveness suitable for multiple choice tests and rating scales.  
Joint entropy $J_{1,2$} is a scalar measure of the cross-product of multinomial 
vectors 1 and 2. Mutual entropy $I_{1,2} = I_1 + I_2 - J_{1,2}$ is a measure 
of the co-dependency of items 1 and 2, and thus the analogue of the negative
log of a squared correlation $R^2$.  this function computes all four types
of entropies for two specificed items.
}
\usage{
entropies(theta, m, n, U, noption)
}
\arguments{
  \item{theta}{A vector of length N containing score index values for each
  test taker.}
  \item{m}{The index of the first choice.}
  \item{n}{The index of the second choice.}
  \item{U}{The data matrix containing the indices of choisen options for
  each test taker.}
  \item{noption}{A vector containing the number of options for all items.}
}
\value{
A named list object containing objects produced from analyzing the simulations,
one set for each simulation:
  \describe{
    \item{I_m:}{The entropy of item m.}
    \item{I_n:}{The entropy of item n.}
    \item{J_nm:}{The joint entropy of items m and n.}
    \item{I_nm:}{The mutual entropy of items m and n.}
  }
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Entropy.plot}}
}
\examples{
#  Load needed objects
U       <- Quantshort_dataList$U
theta   <- Quantshort_parList$theta
noption <- matrix(5,24,1)
#  compute mutual entropies for all pairs of the first 6 items
Mvec    <- 1:6
Mlen    <- length(Mvec)
Hmutual <- matrix(0,Mlen,Mlen)
for (i1 in 1:Mlen) {
  for (i2 in 1:i1) {
    Result <- entropies(theta, Mvec[i1], Mvec[i2], U, noption)
    Hmutual[i1,i2] = Result$Hmutual
    Hmutual[i2,i1] = Result$Hmutual
  }
}
print("Matrix of mutual entries (off-digagonal) and self-entropies (diagonal)")
print(round(Hmutual,3))
}