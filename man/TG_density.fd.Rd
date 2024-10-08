\name{TG_density.fd}
\alias{TG_density.fd}
\title{
  Compute a Probability Density Function
}
\description{
  Like the regular S-PLUS function \code{density}, this function
  computes a probability density function for a sample of values of a
  random variable.  However, in this case the density function is
  defined by a functional parameter object \code{logdensfdPar} along with a
  normalizing constant \code{C}.

  The density function $p(indexdens)$ has the form 
  \code{p(indexdens) = C exp[W(indexdens)]}
  where function $W(indexdens)$ is defined by the functional data object
  \code{logdensfdPar}.
}
\usage{
\method{TG_density}{fd}(indexdens, logdensfd, conv=0.0001, iterlim=20,
           active=1:nbasis, dbglev=0)
}
\arguments{
\item{indexdens}{
a set observations, which may be one of two forms:
\enumerate{
	\item a vector of observations $indexdens_i$
	\item a two-column matrix, with the observations $indexdens_i$ in the
	first column, and frequencies $f_i$ in the second.
}
The first option corresponds to all $f_i = 1$.
}
\item{logdensfd}{
a functional data object specifying the initial
value, basis object, roughness penalty and smoothing
parameter defining function $W(t).$
}
\item{conv}{
a positive constant defining the convergence criterion.
}
\item{iterlim}{
the maximum number of iterations allowed.
}
\item{active}{
a logical vector of length equal to the number of coefficients
defining \code{Wfdobj}. If an entry is TRUE, the corresponding
coefficient is estimated, and if FALSE, it is held at the value defining the
argument \code{Wfdobj}.  Normally the first coefficient is set to 0
and not estimated, since it is assumed that $W(0) = 0$.
}
\item{dbglev}{
either 0, 1, or 2.  This controls the amount information printed out on
each iteration, with 0 implying no output, 1 intermediate output level,
and 2 full output.  If levels 1 and 2 are used, it is helpful to
turn off the output buffering option in S-PLUS.
}
}
\value{
  a named list of length 4 containing:
  \item{Wfdobj}{
    a functional data object defining function $W(indexdens)$ that that
    optimizes the fit to the data of the monotone function that it defines.
  }
  \item{C}{
    the normalizing constant.
  }
  \item{Flist}{
    a named list containing three results for the final converged solution:
    (1)
    \bold{f}: the optimal function value being minimized,
    (2)
    \bold{grad}: the gradient vector at the optimal solution,   and
    (3)
    \bold{norm}: the norm of the gradient vector at the optimal solution.
  }
  \item{iternum}{
    the number of iterations.
  }
  \item{iterhist}{
    a \code{iternum+1} by 5 matrix containing the iteration
    history.
  }
}
\details{
  The goal of the function is provide a smooth density function
  estimate that approaches some target density by an amount that is
  controlled by the linear differential operator \code{Lfdobj} and
  the penalty parameter. For example, if the second derivative of
  $W(t)$ is penalized heavily, this will force the function to
  approach a straight line, which in turn will force the density function
  itself to be nearly normal or Gaussian.  Similarly, to each textbook
  density function there corresponds a $W(t)$, and to each of these
  in turn their corresponds a linear differential operator that will, when
  apply to $W(t)$, produce zero as a result.
  To plot the density function or to evaluate it, evaluate \code{Wfdobj},
  exponentiate the resulting vector, and then divide by the normalizing
  constant \code{C}.
}
\references{
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Full information optimal scoring. 
  Journal of Educational and Behavioral Statistics, 45, 297-315.
  
  Ramsay, J. O., Li J. and Wiberg, M. (2020) Better rating scale scores with 
  information-based psychometrics.  Psych, 2, 347-360.
 }
\author{Juan Li and James Ramsay}
