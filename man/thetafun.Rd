\name{thetafun}
\alias{thetafun}
\title{Compute optimal scores}
\usage{
  thetafun(theta, WfdList, dataList, itermax = 20, crit = 0.001, 
           itdisp = FALSE)
}
\arguments{
  \item{theta}{A vector of size \code{N} containing initial values for score 
  indices in the interval [0,100].}
  \item{WfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  
    Its length is equal to the number of items in the test or questions in the 
    scale.  Each member of \code{WfdList} is a named list containing information 
    computed during the analysis.
    The named list in each member contains these objects:
    \describe{
      \item{Wfd:}{A functional data object containing the \code{M} 
      surprisal curves for a question.}
      \item{M:}{The number of options.}
      \item{Pbin:}{A matrix containing proportions at each bin.}
      \item{Wbin:}{A matrix containing surprisal values at each bin.}
      \item{Pmatfine:}{A matrix of probabilities over a fine mesh.}
      \item{Wmatfine:}{A matrix of surprisal values over a fine mesh.}
      \item{DWmatfine:}{A matrix of the values of the first derivative of 
      surprisal curves over fine mesh.}
      \item{D2Wmatfine:}{A matrix of the values of the second derivative of 
      surprisal curves over fine mesh.}
      \item{chartList:}{A list vector of length \code{n}.  Each list contains 
      the \code{K} by \code{M-1} matrix of initial values for the surprisal 
      curves for the corresponding item.}
    }
  }
  \item{dataList}{A list that contains the objects needed to analyse the test or rating scale
  	with the following fields:
    \describe{
      \item{U:}{A matrix of response data with N rows and n columns where
             N is the number of examinees or respondents and
             n is the number of items.
             Entries in the matrices are the indices of the options chosen.
             Column i of U is expected to contain only the integers 1,...,noption.}  
      \item{optList:}{A list vector containing the numerical score values assigned to the
       options for this question.}
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is a numeric vector 
       	of length n containing the indices the right answers.  Otherwise, it is NULL.}
      \item{chartList:}{A list vector length n each member of which is a matrix of 
      	parameters values defining a set of M surprisal curves.  The matrix has K rows
      	and M - 1 columns where K is the number of basis functions defining the curves and
      	M is the number curves.}
      \item{WfdPar:}{An fdPar object for the defining the surprisal curves.}
      \item{noption:}{A numeric vector of length n containing the numbers of options for 
      	each item.}
      \item{nbin:}{The number of bins for binning the data.}
      \item{scrrng:}{A vector of length 2 containing the limits of observed sum scores.}
      \item{scrfine:}{A fine mesh of test score values for plotting.}
      \item{scrvec:}{A vector of length N containing the examinee or respondent sum scores.}
      \item{itemvec:}{A vector of length n containing the question or item sum scores.}
      \item{percntrnk:}{A vector length N containing the sum score percentile ranks.}
      \item{thetaQnt:}{A numeric vector of length 2*nbin + 1 containing the bin boundaries
      	alternating with the bin centers. These are initially defined as 
      	seq(0,100,len=2*nbin+1).}
      \item{Wdim:}{The total dimension of the surprisal scores.}
      \item{PcntMarkers:}{The marker percentages for plotting: 5, 25, 50, 75 and 95.}
    }
  }
  \item{itermax}{Maximum number of iterations for computing the optimal
  theta values.  Default is 20.}
  \item{crit}{Criterion for convergence of optimization.  Default is 1e-8.}
  \item{itdisp}{If TRUE, results are displayed for each iteration.}
}
\description{
The percentile score index values are estimated for each person.  The estimates
minimize the negative log likelihoods, which are a type of surprisal.  The
main optimization method is a safe-guarded Newton-Raphson method.  

For any iteration the method uses only those scores that are within the interior  
of the interval [0,100] or at a boundary with a first derivative that would
take a step into the interior, and have second derivative values exceeding the 
value of argument \code{crit}.  Consequently the number of values being 
optimized decrease on each iteration, and iterations cease when either 
all values meet the convergence criterion or are optimized on a 
boundary, or when the number of iterations reaches \code{itermax}.
At that point, if there are any interior scores still associated with 
either non-positive second derivatives or values that exceed 
\code{crit}, the minimizing value along a fine mesh is used.

If \code{itdisp} is positive, the number of values to be estimated
are printed for each iteration.
}
\value{
A named list with these members:
  \describe{
    \item{theta_out:}{A vector of optimized score index value.} 
    \item{Hval:}{The negative log likelihood criterion.}
    \item{DHval:}{The first derivative of the negative likelihood.} 
    \item{D2Hval:}{The second derivative of the negative likelihood.} 
    \item{iter:}{The number iterations used.}
  }
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.

http://www.testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{theta.distn}},
  \code{\link{Hfun}},
  \code{\link{DHfun}},
  \code{\link{theta2arclen}},
  \code{\link{testscore}},
  \code{\link{scoreDensity}}
}
\examples{
\donttest{
#  Example 1.  Compute score index values theta for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
theta_in  <- Quant_dataList$percntrnk
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(2,1))
hist(theta_in,11)
WfdList   <- Quant_parListi$WfdList
thetafunResult <- thetafun(theta_in, WfdList, Quant_dataList)
theta_out      <- thetafunResult$theta_out
hist(theta_out,11)
par(oldpar)
#  Example 2.  Compute score index values theta for the 
#  Symptom Distress Scale with 13 items and 473 respondents
#  Proceed as above changing "Quant" to "SDS".
}
}
