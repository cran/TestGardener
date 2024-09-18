\name{Sbinsmth.init}
\alias{Sbinsmth.init}

\title{
Initialize surprisal smoothing of choice data.
}
\description{
This version of Sbinsmth.init() uses direct least squares smoothing of the
surprisal values at bin centers to generate dependent variables for
a model for the vectorized K by M-1 parameter matrix Bmat.  The estimates
of the surprisal curves are approximated using functions in the fda package.}
\usage{
  Sbinsmth.init(percntrnk, nbin, Sbasis, grbgvec, noption, chcemat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{percntrnk}{Percent rank values of sum score values, usually after
  jittering}
  \item{nbin}{The number of bins used to bin the choice data.}
  \item{Sbasis}{A bspline functional basis object for surprisal smoothing.}
  \item{grbgvec}{A logical vector of length n indicating whether or not
  the choice data for an item is added for missing of illigetimate choices.}
  \item{noption}{An integer vector indicating the number of options for
  each item, not including a possible added garbage option.}
  \item{chcemat}{An N by n matrix with each row containing the indices of 
  the options chosen by a person.}
}
\value{A list vector of length n, each element being a list vector containing
objects necessary for surprisal smoothing.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{ICC_plot}},
  \code{\link{Sbinsmth}}
}
