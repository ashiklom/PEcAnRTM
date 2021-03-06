% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invert.custom.R
\name{invert.custom}
\alias{invert.custom}
\title{Bayesian inversion of a model}
\usage{
invert.custom(
  observed,
  invert.options,
  quiet = FALSE,
  return.resume = FALSE,
  runID = NULL
)
}
\arguments{
\item{observed}{Vector, matrix, or data frame (coerced to matrix) of
observed values. For spectral data, wavelengths are rows and spectra are
columns. Dimensions must align with the output of \code{model}.}

\item{invert.options}{R list object containing inversion settings. See details.}

\item{quiet}{Suppress progress bar and status messages. Default=FALSE}

\item{return.resume}{If \code{TRUE}, return results as list that includes current
Jump distribution (useful for continuing an ongoing run) and acceptance
rate. Default = \code{FALSE}.}

\item{runID}{Run-unique ID. Useful for parallel runs. Default=NULL}
}
\description{
Performs an inversion of an arbitrary model using a modified
Metropolis Hastings algorithm with block sampling. This may be slightly
slower than the implementation in Fortran, but is much more customizable, as
the model can be any R function.
}
\details{
\code{inversion.options} contains the following:
\itemize{
\item \code{inits} -- Vector of initial values of model parameters to be inverted.
\item \code{ngibbs} -- Number of MCMC iterations
\item \code{prior.function} -- Function for use as prior.
Should take a vector of parameters as input and return a single value -- the
sum of their log-densities -- as output.
\item \code{param.mins} -- Vector of minimum values for inversion parameters
\item \code{param.maxs} -- Vector of minimum values for inversion parameters
\item \code{model} -- The model to be inverted.
This should be an R function that takes \code{params} and \code{runID} as input and
returns one column of \code{observed} (nrows should be the same).
Constants should be implicitly included here.
\item \code{adapt} -- Number of steps for adapting covariance matrix (i.e. adapt
every 'n' steps). Default=100
\item \code{adj_min} -- Minimum threshold for rescaling Jump standard deviation.
Default = 0.1.
\item \code{target} -- Target acceptance rate. Default=0.234, based on recommendation
for multivariate block sampling in Haario et al. 2001
\item \code{do.lsq} -- Perform least squares optimization first (see \code{invert.lsq}),
and use outputs to initialize Metropolis Hastings.
This may improve mixing time, but risks getting caught in a local minimum.
Default=FALSE
\item \code{catch_error} -- If \code{TRUE} (default), wrap model in \code{tryCatch} to prevent sampling termination on model execution error.
}
}
\references{
\itemize{
\item Haario, Heikki; Saksman, Eero; Tamminen, Johanna.  An adaptive Metropolis
algorithm. Bernoulli 7 (2001), no. 2, 223--242.
http://projecteuclid.org/euclid.bj/1080222083.
}
}
