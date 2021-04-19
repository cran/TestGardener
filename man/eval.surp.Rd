\name{eval.surp}
\alias{eval.surp}
\title{
  Evaluate multinomial surprisal curves at a vector of score index values.
}
\description{
  Surprisal is -log(probability) where the logarithm is equal to the length 
  \code{M} of the multinomial observation vector.  The surprisal curves for each question are estimated
  by fitting the surprisal values of binned data using curves whose values are within the
  M-1 dimensional surprisal subspace that is within the space of non-negative M-dimensional vectors.
  This function computes the values of the \code{M} surprisal curves at a vector of argument values.
}
\usage{
  eval.surp(evalarg, Wfdobj, nderiv=0)
}
\arguments{
  \item{evalarg}{
    Argument value array of length N, where N is the number of observed
    curve values for each curve.  It is assumed that that these argument
    values are common to all observed curves.  If this is not the case,
    you will need to run this function inside one or more loops,
    smoothing each curve separately.
  }
  \item{Wfdobj}{
    A functional data object that defines the surprisal curves.
  }
  \item{nderiv}{
    The order to derivative the curves to be evaluated.  This must be
    0, 1 or 2.
  }
}
\value{
A \code{length(argvals)} by \code{M} matrix of surprisal derivative values.}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with information-based 
psychometrics.  Psych, 2, 347-360.

http://www.testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{smooth.surp}}
}
\examples{
#  see example for smooth.surp()
}
