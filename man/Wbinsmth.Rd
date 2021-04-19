\name{Wbinsmth}
\alias{Wbinsmth}
\title{Estimate the option probability and surprisal curves.}
\usage{
  Wbinsmth(theta, dataList, thetaQnt=seq(0,100,len=2*nbin+1), chartList)
}
\arguments{
  \item{theta}{A vector of length N containing current values of score index 
  percentile values.}  
  \item{dataList}{A list that contains the objects needed to analyse the test or 
  rating scale with the following fields:
    \describe{
      \item{U:}{A matrix of response data with N rows and n columns where
             N is the number of examinees or respondents and
             n is the number of items.
             Entries in the matrices are the indices of the options chosen.
             Column i of U is expected to contain only the integers 
             \code{1,...,noption}.}  
      \item{optList:}{A list vector containing the numerical score values 
      assigned to the options for this question.}
      \item{key:}{If the data are from a test of the multiple choices type
       	where the right answer is scored 1 and the wrong answers 0, this is 
       	a numeric vector of length n containing the indices the right answers.  
       	Otherwise, it is NULL.}
      \item{chartList:}{A list vector length n each member of which is a 
      matrix of parameters values defining a set of M surprisal curves.  
      The matrix has K rows and M - 1 columns where K is the number of basis 
      functions defining the curves and M is the number curves.}
      \item{WfdPar:}{An fdPar object for the defining the surprisal curves.}
      \item{noption:}{A numeric vector of length n containing the numbers of 
      options for each item.}
      \item{nbin:}{The number of bins for binning the data.}
      \item{scrrng:}{A vector of length 2 containing the limits of observed 
      sum scores.}
      \item{scrfine:}{A fine mesh of test score values for plotting.}
      \item{scrvec:}{A vector of length N containing the examinee or 
      respondent sum scores.}
      \item{itemvec:}{A vector of length n containing the question or item 
      sum scores.}
      \item{percntrnk:}{A vector length N containing the sum score 
      percentile ranks.}
      \item{thetaQnt:}{A numeric vector of length 2*nbin + 1 containing 
      the bin boundaries alternating with the bin centers. These are initially 
      defined as \code{seq(0,100,len=2*nbin+1)}.}
      \item{Wdim:}{The total dimension of the surprisal scores.}
      \item{PcntMarkers:}{The marker percentages for plotting: 
      5, 25, 50, 75 and 95.}
    }
  }
  \item{thetaQnt}{A vector of length 2*n+1 containing the sequence of bin boundary and 
  bin centre values.}
  \item{chartList}{A list vector of length n containing for each item an
      	K by M-1 dimension matrix to be used to initialize function smooth.surp().
      	K is the number of basis functions used defining a surprisal curve and M
      	is the number of options.}
}
\value{
  \item{WfdList}{
  A numbered list object produced by a TestGardener analysis of a test.  Its length
  is equal to the number of items in the test or questions in the scale.  
  Each member of \code{WfdList} is a named list containing information computed during the analysis.
  These named lists contain these objects:
    \describe{
      \item{Wfd:}{A functional data object containing the \code{M} surprisal curves
      for a question.}
      \item{M:}{The number of options.}
      \item{Pbin:}{A matrix containing proportions at each bin.}
      \item{Wbin:}{A matrix containing surprisal values at each bin.}
      \item{Pmatfine:}{A matrix of probabilities over a fine mesh.}
      \item{Wmatfine:}{A matrix of surprisal values over a fine mesh.}
      \item{DWmatfine:}{A matrix of the values of the first derivative of 
      surprisal curves over fine mesh.}
      \item{D2Wmatfine:}{A matrix of the values of the second derivative of 
      surprisal curves over fine mesh.}
      \item{chartList:}{A list vector of length \code{n}.  Each list contains the
      \code{K} by \code{M-1} matrix of initial values for the surprisal curves for
      the corresponding item.}
    }
  }
}
\description{
  The function first bins the data in order to achieve rapid estimation of the option 
  surprisal curves. The argument \code{thetaQnt} contains the sequence of bin boundaries 
  separated by the bin centers, so that it is of length \code{2*nbin + 1} where \code{nbin}
  is the number of bins.  These bin values are distributed over the percentile interval
  [0,100] so that the lowest boundary is 0 and highest 100.  Prior to the call to 
  \code{Wbinsmth} these boundaries are computed so that the numbers of values of 
  \code{theta} falling in the bins are roughly equal.  It is important that the
  number of bins be chosen so that the bins contain at least about 25 values.  
  
  After the values of \code{theta} are binned, the proportions that the bins are chosen 
  for each question and each option are computed.  Proportions of zero are given NA 
  values.
  
  The positive proportions are then converted to surprisal values where surprisal
  = -log_M (proportion) where log_M is the logarithm with base M, the number 
  of options associated with a question.  Bins with zero proportions are assigned 
  a surprisal that is appropriately large in the sense of being in the range of
  the larger surprisal values associated with small but positive proportions.
  
  The next step is to fit the surprisal values for each question by a 
  functional data object that is smooth, passes as closely as possible to an
  option's surprisal values, and has values consistent with being a surprisal 
  value.  The function \code{smooth.surp()} is used for this purpose.
  
  Finally the curves and other results for each question are saved in object
  \code{WfdList}, a list vector of length n, and the list vector is returned.
}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics. Psych, 2, 347-360.

http://www.testgardener.azurewebsites.net
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{Wbinsmth.plot}},
  \code{\link{thetafun}}
}
\examples{
#  Example 1.  Display the item probability and surprisal curves for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
#  estimate the initial surprisal curves
theta     <- Quant_dataList$percntrnk
chartList <- Quant_dataList$chartList
thetaQnt  <- Quant_parListi$thetaQnt
Qvec      <- Quant_parListi$Qvec
WfdResult_init <- Wbinsmth(theta, Quant_dataList, thetaQnt, chartList)
WfdList_init   <- WfdResult_init$WfdList
#  plot the curves for the first question
binctr <- Quant_parListi$binctr
Wbinsmth.plot(binctr, Qvec, WfdList_init, Quant_dataList, plotindex=1)
#  estimate the final surprisal curves
theta    <- Quant_parListi$theta
thetaQnt <- Quant_parListi$thetaQnt
Qvec     <- Quant_parListi$Qvec
WfdResult_final <- Wbinsmth(theta, Quant_dataList, thetaQnt, chartList)
WfdList_final   <- WfdResult_final$WfdList
#  plot the curves for the first question
Wbinsmth.plot(binctr, Qvec, WfdList_final, Quant_dataList, plotindex=1)
#  Example 2.  Display the item probability and surprisal curves for the 
# Symptom Distress Scale with 13 items and 473 respondents.
# Proceed as above changing "Quant" to "SDS"
}
