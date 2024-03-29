\name{mu_plot}
\alias{mu_plot}
\title{
  Plot expected test score as a function of score index
}
\description{
   The expected score \code{mu(index)} is a function of the score index index.  A diagonal
   dashed line is displayed to show the linear relationship to the score range interval.
}
\usage{
mu_plot(mufine, scrrng, titlestr)
}
\arguments{
  \item{mufine}{A mesh of 101 equally spaced values of mu as a function of index.}
  \item{scrrng}{A vector of length 2 containing the score range.}
  \item{titlestr}{A string containing the title of the data.}
}
\value{
  A \code{gg} or \code{ggplot} object defining the plot of the expected test score 
  \code{mu} as a function of   the score index \code{index}.  This is displayed by 
  the \code{print} command. 
  The plot is automatically displayed as a side value even if no return object is 
  specified in the calling statement.
}
\references{
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
  Journal of Educational and Behavioral Statistics, 45, 297-315.
  
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
  information-based psychometrics.  Psych, 2, 347-360.
}
\author{Juan Li and James Ramsay}
\seealso{
  \code{\link{scoreDensity},
        \link{mu}
  }
}
