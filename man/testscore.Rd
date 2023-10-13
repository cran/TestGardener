\name{testscore}
\alias{testscore}
\title{Compute the expected test score by substituting probability of choices
for indicator variable 0-1 values.
Binary items assumed coded as two choice items.}
\usage{
  testscore(theta, WfdList, optList)
}
\arguments{
  \item{theta}{Initial values for score indices in the interval [0,100]. 
    A vector of size \code{N}.}
  \item{WfdList}{
    A numbered list object produced by a TestGardener analysis of a test.  Its length
    is equal to the number of items in the test or questions in the scale.  
    Each member of \code{WfdList} is a named list containing information computed 
    during the analysis.}
  \item{optList}{A numbered list of length n.  Each member contains the weights
    assigned to each option for that item or question.}
}
\description{
  Compute the expected test score by substituting probability of choices
  for indicator variable 0-1 values.
  Binary items assumed coded as two choice items.
}
\value{A vector of test score values.}
\author{Juan Li and James Ramsay}
\references{
Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
Journal of Educational and Behavioral Statistics, 45, 297-315.

Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
information-based psychometrics.  Psych, 2, 347-360.

http://testgardener.azurewebsites.net
}
\seealso{
  \code{\link{scoreDensity}}
}
\examples{
#  Example 1.  Compute expected sum score values for the 
#  short SweSAT multiple choice test with 24 items and 1000 examinees
WfdList <- Quant_13B_problem_parList$WfdList
theta   <- Quant_13B_problem_parList$theta
mu      <- testscore(theta, WfdList, Quant_13B_problem_dataList$optList)
par(c(1,1))
hist(mu,11)
#  Example 2.  Compute expected sum score values for the 
#  Symptom Distress Scale with 13 items and 473 respondents
#  Proceed as above changing "Quant" to "SDS".
}
