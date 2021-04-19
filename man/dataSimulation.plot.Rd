\name{dataSimulation.plot}
\alias{dataSimulation.plot}
\title{
Plot root-mean-squared error and bias for \code{theta} and \code{mu(theta)}
}
\description{
This function is called after a call to function \code{dataSimulation} where a 
set of simulated data sets are analyzed.  This function plots the performance of
estimates of score index \code{theta}, expected score \code{mu(theta)}, and the 
observed sum score.

Mean squared error for a function is a curve showing error variation in an estimated
function.  Root-mean-squared error converts this to a standard error curve just as
a standard deviation does for the variance of a mean, and so provides a more 
visually direct sense of how well an estimated curve is defined by the data.

Mean squared error can be split into two orthogonal or independent components:  
(1) sampling variance showing variation around the mean estimate and
(2) squared bias showing systematic variation around the true curve. 
}
\usage{
dataSimulation.plot(simList, Qvec, ttlsz=NULL, axisttl=NULL, axistxt=NULL, lgdlab=NULL)
}
\arguments{
  \item{simList}{A named list containing these objects:
    \describe{
      \item{sumscr:}{A matrix with row dimension \code{ntheta}, the number of population 
      score index values and column dimension \code{nsample}, the number of simulated 
      samples.}
      \item{theta:}{An \code{ntheta} by \code{nsample} of estimated score index values.}
      \item{mu:}{An \code{ntheta} by \code{nsample} of estimated expected score values.}
      \item{al:}{An \code{ntheta} by \code{nsample} of estimated test information curve 
      values.}
      \item{thepop:}{A vector of population score index values.}
      \item{mupop:}{A vector of expected scores computed from the population score 
      index values.}
      \item{alpop:}{A vector of test information values computed from the population score 
      index values.}
      \item{n:}{The number of questions.}
      \item{Qvec:}{The five marker percentile values.}
    }
  }
  \item{Qvec}{Locations on the score index interval [0,100] of the five marker percentages.}
  \item{ttlsz}{Title font size.}
  \item{axisttl}{Axis title font size.}
  \item{axistxt}{Axis text(tick label) font size.}
  \item{lgdlab}{Legend label font size.}
}
\value{
  Two ggplot objects are returned in a named list with these members:
  \describe{
    \item{thetaplot:}{A plot with two panels.   The upper panel plots the 
    root-mean-squared-error of estimation (RMSE) of the score index value \code{theta}
    as a solid line, and also plots the sampling error (square root of the sampling variance)
    as a red dashed line.  The sampling error is necessarily smaller than the RMSE,
    but if it is only a small amount less, then the bias of the estimate is negligible.
    The lower plot shows the bias in the estimate as a function of \code{theta}.}
    \item{muplot:}{A plot with two panels.  The upper panel plots the 
    root-mean-squared-error of estimation (RMSE) of the expected score value \code{mu(theta)}
    as a function of the score index value \code{theta} as a solid line, and also plots
    the RMSE of the observed sum score as a dashed blue line.  The lower plot shows the bias 
    in these estimate as a function of \code{theta}.
    }
  }
  These plots are also displayed automatically as a side effect.
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
  \code{\link{Usimulate},
        \link{scorePerformance},
        \link{dataSimulation}
  }
}
\examples{
\donttest{
#  Use of dataSimulations() presumes that make.data() has already been executed
#  to set up the analysis, Analyze() has been run with the desired number of
#  iterations and its two objects extracted.
#  Example 1:  Ten simulations of multiple choice test, 24 items, 1000 examinees
Quant_simList <- dataSimulation(Quant_dataList, Quant_parListi, nsample=10)
dataSimulation.plot(Quant_simList, Quant_parListi$Qvec)
#  Example 2:  Ten simulations of 13-item Symptom Distress Scale with 473 respondents
SDS_simList <- dataSimulation(SDS_dataList, SDS_parListi, nsample=10)
dataSimulation.plot(SDS_simList, SDS_parListi$Qvec)
}
}

